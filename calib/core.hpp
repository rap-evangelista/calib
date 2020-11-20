#ifndef DEFAULT_SPACE_DIM
    #define DEFAULT_SPACE_DIM 2
#endif

#ifndef CALIB_MODE_HOST_DEVICE
    #define CALIB_MODE_HOST_ONLY
#endif

#include "datastructure.hpp"

#ifndef CALIB_METRIC_SPACE
    #define CALIB_METRIC_SPACE
    #include "geometries\euclidian.hpp"
#endif