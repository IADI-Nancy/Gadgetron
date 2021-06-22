/*
    Karyna ISAIEVA, 2020
*/
#pragma once

#include <Gadget.h>
#include <hoNDArray.h>
#include "gadgetron_iadi_export.h"

#include <mri_core_data.h>
#include <utility>
#include <vector>

namespace Gadgetron{

  class EXPORTGADGETSIADI ReorderSlicesGadget : 
  public Gadget1Of2<IsmrmrdImageArray, ISMRMRD::ImageHeader >
    {
    public:
      GADGET_DECLARE(ReorderSlicesGadget)
      ReorderSlicesGadget();
	
    protected:
      virtual int process(GadgetContainerMessage<IsmrmrdImageArray>* m1);
      virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1);

    private:
        static bool compare_slices(const std::pair<ISMRMRD::ImageHeader*, std::complex<float>*>& a, 
            const std::pair<ISMRMRD::ImageHeader*, std::complex<float>*>& b);
    };
}