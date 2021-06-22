/*
    Karyna ISAIEVA, 2020
*/
#pragma once

#include <gadgetron/Node.h>
#include <gadgetron/hoNDArray.h>
#include <ismrmrd/ismrmrd.h>

namespace Gadgetron
{
    typedef Gadgetron::hoNDArray<unsigned int> hoArrayUInt;
    typedef std::pair<Gadgetron::hoNDArray<unsigned int>, Gadgetron::hoNDArray<unsigned int>> TimeDataPair_t;
    typedef std::pair<size_t, size_t> SliSet_t;

    struct AcquisitionBuffer {
        hoNDArray<std::complex<float>> data;
        std::vector<ISMRMRD::AcquisitionHeader> headers;
    };

    struct PhysData {
        std::vector<unsigned long> timestamps;
        std::vector<double> val;
    }; 

    struct AcquisitionInfo {
        std::vector<unsigned long> timestamps;
        std::vector<size_t> sli;
        std::vector<size_t> set;
    };

    class SaveKspaceData : public Core::ChannelGadget<Core::variant<Core::Acquisition, Core::Waveform>>
    {

        public:
            using Core::ChannelGadget<Core::variant<Core::Acquisition, Core::Waveform>>::ChannelGadget;
            void process(Core::InputChannel<Core::variant<Core::Acquisition, Core::Waveform>>& in,
                Core::OutputChannel& out) override;

        NODE_PROPERTY(GRICS_folder, std::string, "Path to GRICS binary files", "/opt/data/GRICS");
        NODE_PROPERTY(savePhysiological, bool, "Whether to save Siemens physiological data", true);
        // NODE_PROPERTY(MRI_SAEC_offset, double, "Time between sequence start and the first calib line", 0.);

        private:
            const char *folderNames = "%s/Siemens_SingleImage_slice%02d_image%02d";
            // Map by type of data (e.g. 'ECG' or 'RESP'), the values is a pair of timestamps/phys_datas
            std::map<std::string, TimeDataPair_t> process_waveform(const std::vector<Core::Waveform>& waveforms);
            hoArrayUInt upsample_timestamps(const std::vector<ISMRMRD::WaveformHeader>& heads, unsigned long length);
            hoArrayUInt filter_physiological_data(const hoArrayUInt& in, int order, double fcut, double fps) ;
            void save_phys_data(const std::map<std::string, TimeDataPair_t>& phys_data, const std::string& data_type);
            std::string get_folder_name(size_t sli, size_t set);
            void save_acquisition_timestamps(const AcquisitionInfo &acq_info);
            void save_mri_saec_offset(const unsigned long& first_kspace_timestamp);

            unsigned long sequence_start = 0;
            const double dt = 2.5e-3; // Siemens physiological data timestamp tick
            std::vector<std::string> phys_data_names = {"ECG", "PULS", "RESP", "EXT1", "EXT2",
                                                        "ECG", "PULS", "RESP", "EXT1", "EXT2"};
            const double MRI_SAEC_offset = 0.;
    };

    GADGETRON_GADGET_EXPORT(SaveKspaceData)
}

std::vector<double> Butterworth(const std::vector<double>& indata, double deltaTimeinsec, double CutOff);
std::vector<double> interp1( std::vector<double> &x, std::vector<double> &y, std::vector<double> &x_new);
int findNearestNeighbourIndex( double value, std::vector<double> &x );