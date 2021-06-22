#ifndef GADGETRON_IADI_EXPORT_H_
#define GADGETRON_IADI_EXPORT_H_

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_IADI__)
        #define EXPORTGADGETSIADI __declspec(dllexport)
    #else
        #define EXPORTGADGETSIADI __declspec(dllimport)
    #endif
#else
    #define EXPORTGADGETSIADI
#endif

#endif /* GADGETRON_IADI_EXPORT_H_ */
