#
#    Karyna ISAIEVA, 2020
#
import gadgetron
import ismrmrd
import numpy as np
import sys
import struct
sys.path.insert(0, "/opt/code/iadi/h5Wrapper/")
from h5Saec import *
# from scipy.signal import butter, filtfilt
# import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os
from pyftpdlib.authorizers import DummyAuthorizer
from pyftpdlib.handlers import FTPHandler
from pyftpdlib.servers import FTPServer
import requests
# from scipy.optimize import curve_fit

# To debug the Gadgetron reconstruction offline with .h5 files
read_debug = True
path = '/opt/saec/' # Path inside the docker container
data_path = '/opt/data/GRICS/'

#
offset_filename = "/opt/data/GRICS/saec_mri_offset.dat"
acquisition_info_filename = "/opt/data/GRICS/acquisition_info.dat"
saec_filename_file = path + "saec_filename.txt"

#SAEC Server is configured to send record file over FTP when record is finished
url_keep_alive = 'http://192.168.139.173:8814/records/keepAlive'
url_stop_record = 'http://192.168.139.173:8814/record/stop'

def read_acquisition_info():
    acquisition_timestamps =[] 
    slices = []
    sets = []

    with open(acquisition_info_filename, "r") as timestamps_file:
        for line in timestamps_file:
            values =  line.split('\t')
            timestamp = np.uint32(values[0])
            sli = np.uint8(values[1])
            set = np.uint8(values[2])
            acquisition_timestamps.append(timestamp)
            slices.append(sli)
            sets.append(set)
            
    return acquisition_timestamps, slices, sets

def get_mri_timestamps_in_sec(acquisition_timestamps):
    with open(offset_filename, "r") as offset_file:
        saec_mri_offset = np.float64(offset_file.read())
    timestamps_in_sec = []
    for timestamp in acquisition_timestamps:
        time_in_sec = (np.float64(timestamp) - np.float64(acquisition_timestamps[0])) * 2.5e-3 + saec_mri_offset
        timestamps_in_sec.append(time_in_sec)
    return np.asarray(timestamps_in_sec)

def get_respiration_from_saec(filename):
    SAECData = h5Saec.from_file(filename.strip())
    ticksTo1s = SAECData.attributes.ticksTo1s
    for attr, value in SAECData.__dict__.items():
        if 'SAEC_RESP' in attr:
            respiratory_data = value.RESP.datas.values
            timestampsSAEC = value.RESP.timestamp.values
        if 'SAEC_TRIGGER_SIEMENS' in attr:
            sequence_start = value.SeqStart.timestamp.values[-1]
    timestamps_in_sec = []
    for timestamp in timestampsSAEC:
        timestamp_in_sec = (np.float64(timestamp) - np.float64(sequence_start)) / ticksTo1s
        timestamps_in_sec.append(timestamp_in_sec)
    return np.asarray(timestamps_in_sec), respiratory_data.astype(np.float64)

def get_filtered_and_interpolated_resp_data(timestamps, respiratory_data, timestamps_mri):
    # Choose the track containing the respiratory data (the Maglife belt has 2 tracks, and in our case only one was used)
    sigma1 = np.std(respiratory_data[:, 0])
    sigma2 = np.std(respiratory_data[:, 1])
    idx_track = 0 if sigma1 > sigma2 else 1

    # Some pre-processing here

    # Data resampling
    respiratory_data_f = interp1d(timestamps, respiratory_data[:, idx_track])
    respiratory_data_interpolated = respiratory_data_f(timestamps_mri)
    respiratory_data_interpolated = respiratory_data_interpolated.astype(np.float32)
    respiratory_data_interpolated = (respiratory_data_interpolated - np.mean(respiratory_data_interpolated)) / np.std(respiratory_data_interpolated)
    return respiratory_data_interpolated

def reshape_and_save_resp_data(respiratory_data_interpolated, slices, sets, N_SLI, N_SET):
    files = []
    for i_sli in range(N_SLI):
        for i_set in range(N_SET):
            folder = data_path + 'Siemens_SingleImage_slice' + str(i_sli + 1).zfill(2) + '_image' + str(i_set + 1).zfill(2)
            if not os.path.exists(folder):
                os.mkdir(folder)
            files.append(open(folder + '/ModelInputs.dat','wb'))

    for i in range(len(respiratory_data_interpolated)):
        idx = slices[i] * N_SET + sets[i]
        files[idx].write(bytearray(np.float32(respiratory_data_interpolated[i])))

    for file in files:
        file.close()

def ReceiveSAECDataGadget(connection):
    if not read_debug:
        print('Keep SAEC recording URL request')
        requests.post(url_keep_alive)
        # SAEC Server is configured to cancel record if Gadgetron dont tell that we want to keep this record
        
        # Initialize the server during the gadget pre-load
        authorizer = DummyAuthorizer()
        authorizer.add_user("user", "12345", path, perm="elradfmwMT")
        handler = MyHandler
        handler.authorizer = authorizer
        server = FTPServer(("0.0.0.0", 8021), handler)

    # Receive all the data - syncronization, the gadgets will not be executed until it is done
    acquisitions = []
    for acquisition in connection:
        acquisitions.append(acquisition)

    if not read_debug:
        print('Stopping SAEC recording')        
        myobj = {'sequenceName': connection.header.measurementInformation.protocolName}
        requests.post(url_stop_record, data = myobj)
        # Call HttpClientHandler::postRecordStop in SAEC Server

        # Launch the ftp-server and wait until a file is received
        server.serve_forever()

    N_SLI = connection.header.encoding[0].encodingLimits.slice.maximum + 1
    try:
        N_SET = connection.header.encoding[0].encodingLimits.set.maximum + 1
    except:
        N_SET = 1 

    acquisition_timestamps, slices, sets = read_acquisition_info()
    timestamps_mri = get_mri_timestamps_in_sec(acquisition_timestamps)
    with open(saec_filename_file, "r") as filename_file:
        saec_filename = filename_file.read()
    timestamps_saec, respiratory_data_saec = get_respiration_from_saec(path + saec_filename)
    respiratory_data_interpolated = get_filtered_and_interpolated_resp_data(timestamps_saec, respiratory_data_saec, timestamps_mri)
    reshape_and_save_resp_data(respiratory_data_interpolated, slices, sets, N_SLI, N_SET)

    # Send the unchanged MRI data down the chain
    for acquisition in acquisitions:
        connection.send(acquisition)

class MyHandler(FTPHandler):    
    def on_connect(self):
        print("%s:%s connected" % (self.remote_ip, self.remote_port))

    def on_file_received(self, file):
        print('A file is received\n')
        save_filename(file)
        self.server.close_when_done()
        pass

    def on_incomplete_file_sent(self, file):
        # do something when a file is partially sent
        pass

    def on_incomplete_file_received(self, file):
        # remove partially uploaded files
        import os
        os.remove(file)

def save_filename(saec_filename):
    saec_filename = saec_filename.split('/')[-1]
    f = open(path + saec_filename_file, 'w')
    f.write(saec_filename)
    f.close()

if __name__ == '__main__':
    gadgetron.external.listen(18003, ReceiveSAECDataGadget)