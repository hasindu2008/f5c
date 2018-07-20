#ifndef FAST5LITE_H
#define FAST5LITE_H

#include <hdf5/serial/hdf5.h>


//#define DEBUG_FAST5_IO 1
#define RAW_ROOT "/Raw/Reads/"



struct fast5_t{
	float* rawptr;
	hsize_t nsample;


    //  Information for scaling raw data from ADC values to pA (are these duplicates?)
    float digitisation;
    float offset;
    float range;
    float sample_rate;

};

typedef struct fast5_t fast5;


// // From scrappie
// typedef struct {
//     //  Information for scaling raw data from ADC values to pA
//     float digitisation;
//     float offset;
//     float range;
//     float sample_rate;
// } fast5_raw_scaling;





inline hid_t fast5_open(char*  filename)
{
    hid_t hdf5file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    return hdf5file;
}

//
inline void fast5_close(hid_t hdf5_file)
{
    H5Fclose(hdf5_file);
}


// //
// std::string fast5_get_raw_read_group(hid_t hdf5_file)
// {
//     std::string read_name = fast5_get_raw_read_name(hdf5_file);
//     if(read_name != "") {
//         return std::string(RAW_ROOT) + read_name;
//     } else {
//         return "";
//     }
// }

// //
// std::string fast5_get_raw_read_name(hid_t hdf5_file)
// {
//     // This code is From scrappie's fast5_interface

//     // retrieve the size of the read name
//     ssize_t size =
//         H5Lget_name_by_idx(hdf5_file, RAW_ROOT, H5_INDEX_NAME, H5_ITER_INC, 0, NULL,
//                            0, H5P_DEFAULT);

//     if (size < 0) {
//         return "";
//     }

//     // copy the read name out of the fast5
//     char* name = (char*)calloc(1 + size, sizeof(char));
//     H5Lget_name_by_idx(hdf5_file, RAW_ROOT, H5_INDEX_NAME, H5_ITER_INC, 0, name,
//                        1 + size, H5P_DEFAULT);

//     // cleanup
//     std::string out(name);
//     free(name);
//     return out;
// }

// //
// std::string fast5_get_read_id(hid_t hdf5_file)
// {
//     int ret;
//     hid_t read_name_attribute, raw_group, attribute_type;
//     size_t storage_size = 0;
//     char* read_name_str = NULL;

//     std::string out = "";
    
//     // Get the path to the raw read group
//     std::string raw_read_group = fast5_get_raw_read_group(hdf5_file);



//     if(raw_read_group == "") {
//         return out;
//     }

//     return fast5_get_fixed_string_attribute(hdf5_file, raw_read_group, "read_id");
// }

//inline correct?
// from scrappie
inline float fast5_read_float_attribute(hid_t group, const char *attribute) {
    float val = NAN;
    if (group < 0) {
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "Invalid group passed to %s:%d.", __FILE__, __LINE__);
#endif
        return val;
    }

    hid_t attr = H5Aopen(group, attribute, H5P_DEFAULT);
    if (attr < 0) {
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "Failed to open attribute '%s' for reading.", attribute);
#endif
        return val;
    }

    H5Aread(attr, H5T_NATIVE_FLOAT, &val);
    H5Aclose(attr);

    return val;
}




//
inline int fast5_read(hid_t hdf5_file, fast5* fast5file)
{

    fast5file->rawptr = NULL;
    hid_t space;
    hsize_t nsample;
    herr_t status;


    // retrieve the size of the read name
    ssize_t size =H5Lget_name_by_idx(hdf5_file, RAW_ROOT, H5_INDEX_NAME, H5_ITER_INC, 0, NULL,0, H5P_DEFAULT);

    if (size < 0) {
        fprintf(stderr,"How comw?\n");
    }

    // copy the read name out of the fast5
    char* read_name = (char*)calloc(1 + size, sizeof(char));
    H5Lget_name_by_idx(hdf5_file, RAW_ROOT, H5_INDEX_NAME, H5_ITER_INC, 0, read_name,
                       1 + size, H5P_DEFAULT);

    //not the most efficient and safest, but being a bit lazy for the moment
    char *signal_path=(char *)malloc((size+1+strlen(RAW_ROOT)+strlen("/Signal"))*sizeof(char));
    sprintf(signal_path,"%s%s%s",RAW_ROOT,read_name,"/Signal");	
    printf("Signal path : %s\n",signal_path);
    free(read_name);

    hid_t dset = H5Dopen(hdf5_file, signal_path, H5P_DEFAULT);
    if (dset < 0) {
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "Failed to open dataset '%s' to read raw signal from.\n", signal_path);
#endif
        goto cleanup2;
    }

    space = H5Dget_space(dset);
    if (space < 0) {
        fprintf(stderr, "Failed to create copy of dataspace for raw signal %s.\n", signal_path);
        goto cleanup3;
    }

    H5Sget_simple_extent_dims(space, &nsample, NULL);
    fast5file->nsample=nsample;
    fast5file->rawptr = (float*)calloc(nsample, sizeof(float));
    status = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, fast5file->rawptr);

    if (status < 0) {
        free(fast5file->rawptr);
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "Failed to read raw data from dataset %s.\n", signal_path);
#endif
        goto cleanup4;
    }


 cleanup4:
    H5Sclose(space);
 cleanup3:
    H5Dclose(dset);
 cleanup2:


    //get channel parameters
    const char *scaling_path = "/UniqueGlobalKey/channel_id";

    hid_t scaling_group = H5Gopen(hdf5_file, scaling_path, H5P_DEFAULT);
    if (scaling_group < 0) {
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "Failed to group %s.", scaling_path);
#endif
        return 0;
    }

    fast5file->digitisation = fast5_read_float_attribute(scaling_group, "digitisation");
    fast5file->offset = fast5_read_float_attribute(scaling_group, "offset");
    fast5file->range = fast5_read_float_attribute(scaling_group, "range");
    fast5file->sample_rate = fast5_read_float_attribute(scaling_group, "sampling_rate");

    H5Gclose(scaling_group);


    return 1;

}


#endif