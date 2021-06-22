/*
    Karyna ISAIEVA, 2020
*/
#include "ReorderSlicesGadget.h"
#include <algorithm>

namespace Gadgetron{

ReorderSlicesGadget::ReorderSlicesGadget()
{

}

int ReorderSlicesGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1)
{
    if (this->next()->putq(m1) < 0)
    {
        m1->release();
        return GADGET_FAIL;
    }

    return GADGET_OK;
}

int ReorderSlicesGadget::process( GadgetContainerMessage<IsmrmrdImageArray>* m1)
{
    
    //Grab a reference to the buffer containing the imaging data
    IsmrmrdImageArray& imagearr = *m1->getObjectPtr();

    //7D, fixed order [X, Y, Z, CHA, N, S, LOC]
    uint16_t X = imagearr.data_.get_size(0);
    uint16_t Y = imagearr.data_.get_size(1);
    uint16_t Z = imagearr.data_.get_size(2);
    uint16_t CHA = imagearr.data_.get_size(3);
    uint16_t N = imagearr.data_.get_size(4);
    uint16_t S = imagearr.data_.get_size(5);
    uint16_t LOC = imagearr.data_.get_size(6);

    //Each image will be [X,Y,Z,CHA] big
    std::vector<size_t> img_dims(4);
    img_dims[0] = X;
    img_dims[1] = Y;
    img_dims[2] = Z;
    img_dims[3] = CHA;

    std::vector<std::pair<ISMRMRD::ImageHeader*, std::complex<float>* > > imagearr_vec;

    //Loop over N, S and LOC
    for (uint16_t loc=0; loc < LOC; loc++) {
        for (uint16_t s=0; s < S; s++) {                
            for (uint16_t n=0; n < N; n++) {
                imagearr_vec.push_back(std::make_pair(&imagearr.headers_(n,s,loc), &imagearr.data_(0,0,0,0,n,s,loc)));
            }
        }
    }
    std::sort(imagearr_vec.begin(), imagearr_vec.end(), ReorderSlicesGadget::compare_slices);

    // Prepare the output package
    GadgetContainerMessage<IsmrmrdImageArray>* cm1 = 
                        new GadgetContainerMessage<IsmrmrdImageArray>();
    auto imarray = cm1->getObjectPtr();
    imarray->data_.create({X, Y, Z, CHA, N, S, LOC});
    imarray->headers_.create({N, S, LOC});
    for (uint16_t loc=0; loc < LOC; loc++) {
        for (uint16_t s=0; s < S; s++) {                
            for (uint16_t n=0; n < N; n++) {
                imagearr_vec[loc].first->slice = loc;
                imagearr_vec[loc].first->image_index = loc + 1;
                memcpy(&imarray->headers_(n,s,loc), imagearr_vec[loc].first, sizeof(ISMRMRD::ImageHeader));
                memcpy(&imarray->data_(0,0,0,0,n,s,loc), imagearr_vec[loc].second, X*Y*Z*CHA*sizeof(std::complex<float>));
            }
        }
    }  

    if (this->next()->putq(cm1) < 0) {
        m1->release();
        return GADGET_FAIL;
    } 
    
    m1->release();
    return GADGET_OK;  

}

bool ReorderSlicesGadget::compare_slices(const std::pair<ISMRMRD::ImageHeader*, std::complex<float>*>& a, 
            const std::pair<ISMRMRD::ImageHeader*, std::complex<float>*>& b)
{
    // Works only for axial slices in HFS position
    // Todo : other slice and patient positions
    double slice_a = a.first->position[2];
    double slice_b = b.first->position[2];
    return slice_a > slice_b;
}

GADGET_FACTORY_DECLARE(ReorderSlicesGadget)
}