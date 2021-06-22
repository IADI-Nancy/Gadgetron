PROJECT (gadgetron_iadi LANGUAGES CXX C)

CMAKE_MINIMUM_REQUIRED (VERSION 3.10...3.12 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

### General definitions ###
SET(GADGETRON_INCLUDE_DIR /usr/local/include/gadgetron/)
SET(GADGETRON_LIB_DIR /usr/local/lib/)
SET(GADGETRON_HOME /usr/local/)
SET(GADGETRON_PYTHON /usr/local/share/gadgetron/python/)
SET(GADGETRON_CONFIG_PATH /usr/local/share/gadgetron/config/)
SET(GADGETRON_CODE_DIR /opt/code/)

SET(CMAKE_INSTALL_PREFIX ${GADGETRON_HOME})

FIND_PACKAGE(Boost COMPONENTS thread REQUIRED)
MESSAGE("Boost found: " ${Boost_FOUND})


INCLUDE_DIRECTORIES(
    ${ISMRMRD_INCLUDE_DIR}
    ${HDF5_C_INCLUDE_DIR}
    ${GADGETRON_INCLUDE_DIR}
    # /usr/include/opencv4/
    ${CURL_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
)

### Project properties ###
SET(CMAKE_BUILD_TYPE, "Release")

LINK_DIRECTORIES(${GADGETRON_LIB_DIR})

SET (gadgetron_iadi_header_files 
    gadgetron_iadi_export.h
    SaveKspaceData.h
    ReorderSlicesGadget.h
    )

SET (gadgetron_iadi_src_files 
    SaveKspaceData.cpp
    ReorderSlicesGadget.cpp
    )

FILE (GLOB config_iadi_files
    "config/*.xml")

SET (python_iadi_files
    ReceiveSAECDataGadget.py
    )

ADD_LIBRARY (gadgetron_iadi SHARED 
    ${gadgetron_iadi_header_files} 
    ${gadgetron_iadi_src_files}
    ${config_iadi_files} )

SET_TARGET_PROPERTIES (gadgetron_iadi PROPERTIES
                        VERSION ${GADGETRON_VERSION_STRING}
                        SOVERSION ${GADGETRON_SOVERSION})

TARGET_LINK_LIBRARIES (gadgetron_iadi
                        ${Boost_LIBRARIES}
)

INSTALL (FILES ${gadgetron_iadi_header_files} DESTINATION ${GADGETRON_INCLUDE_DIR})
INSTALL (FILES ${config_iadi_files} DESTINATION ${GADGETRON_CONFIG_PATH})
INSTALL (FILES ${python_iadi_files} DESTINATION ${GADGETRON_PYTHON})

INSTALL (TARGETS gadgetron_iadi
        DESTINATION lib
)