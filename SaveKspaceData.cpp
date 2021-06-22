/*
    Karyna ISAIEVA, 2020
*/
#include "SaveKspaceData.h"
#include <map>
#include "mri_core_data.h"
#include <fstream>
#include <hoNDArray_utils.h>
#include <hoNDImage_util.h>
#include <math.h>
#include <boost/filesystem.hpp>

using namespace Gadgetron;

bool is_noise(Core::Acquisition& acq) {
    return std::get<ISMRMRD::AcquisitionHeader>(acq).isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT);
}

std::string SaveKspaceData::get_folder_name(size_t sli, size_t set) {
    char buffer[100];
    std::sprintf(buffer, folderNames, GRICS_folder.c_str(),
    int(sli + 1), int(set + 1));
    return std::string(buffer);
}

void SaveKspaceData::process(Core::InputChannel<Core::variant<Core::Acquisition, Core::Waveform>>& in,
    Core::OutputChannel& out)
{
    if (!boost::filesystem::exists(GRICS_folder))
        boost::filesystem::create_directory(GRICS_folder);
    
    auto waveforms = std::vector<Core::Waveform>{};
    std::map<SliSet_t, std::ofstream*> kspace_files;
    std::map<SliSet_t, std::ofstream*> label_files;
    std::map<SliSet_t, std::vector<unsigned int>> acquisition_time_stamps;

    //PhysData phys_data;
    AcquisitionInfo acq_info;
    acq_info.timestamps = std::vector<unsigned long>();
    acq_info.sli = std::vector<size_t>();
    acq_info.set = std::vector<size_t>();

    for (auto message : in)
    {
        // Store the waveform messages for the further post-processing
        if (Core::holds_alternative<Core::Waveform>(message)) {
            auto waveform = Core::get<Core::Waveform>(message);
            waveforms.emplace_back(std::move(waveform));
            out.push(message);
            continue;
        }
        auto& acq = Core::get<Core::Acquisition>(message);
        auto head = std::get<ISMRMRD::AcquisitionHeader>(acq);
        auto data = std::get<Gadgetron::hoNDArray<std::complex<float>>>(acq);   

        if (sequence_start == 0)
            sequence_start = head.acquisition_time_stamp;   

        if (is_noise(acq)) {
            out.push(message);
            continue;
        }        

        size_t e1 = head.idx.kspace_encode_step_1;
        size_t sli = head.idx.slice;
        size_t set = head.idx.set;
        SliSet_t sli_set = {sli, set};

        // Check if the file for the pair of slice/set exists in the map
        bool file_not_exists = (kspace_files.find(sli_set) == kspace_files.end());
        // If does not exist, create and open for both Kspace data and Kspace labels
        if (file_not_exists) {
            auto dir_name = get_folder_name(sli, set);
            if (!boost::filesystem::exists(dir_name))
                boost::filesystem::create_directory(dir_name);
            kspace_files[sli_set] = new std::ofstream(dir_name + "/KspaceData.dat", std::ios::out | std::ios::binary);
            label_files[sli_set] = new std::ofstream(dir_name + "/KspaceLabels.dat", std::ios::out | std::ios::binary);
        }
        // Write the kspace line into the file
        std::ofstream* kspace_file = kspace_files.at(sli_set);
        for (size_t c = 0; c < data.get_size(1); c++) {
            for (size_t e0 = 0; e0 < data.get_size(0); e0++) {
                double real_val = double(data(e0, c).real());
                double imag_val = double(data(e0, c).imag());
                kspace_file->write((char *) &real_val, sizeof(double));
                kspace_file->write((char *) &imag_val, sizeof(double));
            }
        }
        // Write index of the kspace line to the file
        std::ofstream* label_file = label_files.at(sli_set);
        int e1_int = e1;
        label_file->write((char *) &e1_int, sizeof(int));
        acquisition_time_stamps[sli_set].push_back(head.acquisition_time_stamp);

        if (acq_info.timestamps.size() == 0) {
          std::cout << "First timestamp : " << head.acquisition_time_stamp << std::endl;
          unsigned long temp = head.acquisition_time_stamp;
          std::cout << "First timestamp : " << temp << std::endl;
          save_mri_saec_offset(temp);
        }
        acq_info.timestamps.push_back(head.acquisition_time_stamp);
        acq_info.sli.push_back(sli);
        acq_info.set.push_back(set);

        out.push(message);
    }
    // Close and delete the files
    for (auto it = kspace_files.begin(); it != kspace_files.end(); it++) {
        it->second->close();
        delete it->second;
    }
    for (auto it = label_files.begin(); it != label_files.end(); it++) {
        it->second->close();
        delete it->second;
    }

    if (savePhysiological) {
        // Rearrange and filter the physiological data
        auto phys_data = process_waveform(waveforms);
        // Save the selected data for the GRICS reconstruction
        save_phys_data(phys_data, "RESP");
    }
    save_acquisition_timestamps(acq_info);
}

void SaveKspaceData::save_acquisition_timestamps(const AcquisitionInfo &acq_info)
{
    std::ofstream file(GRICS_folder + "/acquisition_info.dat");
    for (int i = 0; i < acq_info.timestamps.size(); i++) {
        size_t sli = acq_info.sli[i];
        size_t set = acq_info.set[i];
        auto timestamp = acq_info.timestamps[i]; 
        file << timestamp << "\t" << sli << "\t" << set << std::endl;
    }
    file.close();
}

void SaveKspaceData::save_mri_saec_offset(const unsigned long& first_kspace_timestamp)
{
    std::cout << "Offset info: " << sequence_start << "\t" << first_kspace_timestamp << "\t" << MRI_SAEC_offset << std::endl;
    unsigned long timestamp_dif = first_kspace_timestamp - sequence_start;
    double offset = double(timestamp_dif) * dt + MRI_SAEC_offset;
    std::ofstream file(GRICS_folder + "/saec_mri_offset.dat");
    file << offset << std::endl;
    file.close();
}


std::map<std::string, TimeDataPair_t> SaveKspaceData::process_waveform(const std::vector<Core::Waveform>& waveforms)
{
    std::cout << "Start processing phys data" << std::endl;
    std::map<std::string, std::pair<std::vector<ISMRMRD::WaveformHeader>, std::vector<hoArrayUInt>>> datas;
    for (auto it = waveforms.begin(); it != waveforms.end(); it++) {
        auto head = std::get<ISMRMRD::WaveformHeader>(*it);
        auto info = phys_data_names[head.waveform_id];
        auto data = std::get<hoArrayUInt>(*it);
        datas[info].first.push_back(head);
        datas[info].second.push_back(data);
    }

    std::map<std::string, TimeDataPair_t> phys_data;
    for (auto it = datas.begin(); it != datas.end(); it++) {
        phys_data[it->first].second = Gadgetron::concat_along_dimension(it->second.second, 0);
        // Some pre-processing here
        const size_t N = phys_data[it->first].second.get_size(0);
        phys_data[it->first].first = upsample_timestamps(it->second.first, N);
    }
    std::cout << "Phys data is processed" << std::endl;
    return phys_data;
}

void SaveKspaceData::save_phys_data(const std::map<std::string, TimeDataPair_t>& phys_data, const std::string& data_type)
{
    //std::ofstream file(dir_name + "/ModelInputs.dat", std::ios::out | std::ios::binary);
    // Todo : separation on different slice and set
    std::ofstream file(GRICS_folder + "/phys_data.dat");
    std::cout << "Saving phys data" << std::endl;
    auto data = phys_data.at(data_type);
    for (int i = 0; i < data.first.get_size(0); i++) {
        file << data.first(i);
        for (int j = 0; j < data.second.get_size(1); j++)
            file << "\t" << data.first(i, j);
        file << std::endl; 
    }
    file.close();
}

hoArrayUInt SaveKspaceData::upsample_timestamps(const std::vector<ISMRMRD::WaveformHeader>& heads,
    unsigned long length)
{    
    long ts_dif = heads[1].time_stamp - heads[0].time_stamp;
    auto ts_step = ts_dif / heads[0].number_of_samples;
    auto waveform_time_stamps = hoArrayUInt({length});
    size_t i_new = 0;
    for (size_t i_old = 0; i_old < heads.size(); i_old++) {
        waveform_time_stamps(i_new) = heads[i_old].time_stamp;
        for (size_t shift = 1; shift < heads[i_old].number_of_samples; shift++) {
            i_new++;
            waveform_time_stamps(i_new) = waveform_time_stamps(i_new - 1) + ts_step;
        }
        i_new++;
    }
    return waveform_time_stamps;
}