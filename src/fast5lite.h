/* @file fast5lite.h
**
** lightweight fast5 reading
** adapted from Nanopolish: nanopolish_fast5_io.cpp originally authored by Jared Simpson
** Code was adapted by Hasindu Gamaarachchi
** @@
******************************************************************************/

#ifndef FAST5LITE_H
#define FAST5LITE_H

#ifndef DISABLE_HDF5

#ifndef HAVE_CONFIG_H
#    define HAVE_CONFIG_H
#    include "config.h"
#endif

#ifdef HAVE_HDF5_SERIAL_HDF5_H
#    include <hdf5/serial/hdf5.h>
#endif

#ifdef HAVE_HDF5_H
#    include <hdf5.h>
#endif

#ifdef HAVE_HDF5_HDF5_H
#    include <hdf5/hdf5.h>
#endif

#ifdef HAVE___HDF5_INCLUDE_HDF5_H
#    include <hdf5.h>
#endif

#else
    typedef int hid_t;
#endif

#include "error.h"
#include <string>

typedef struct{
    hid_t hdf5_file;
    bool is_multi_fast5;
}  fast5_file_t;


#ifndef DISABLE_HDF5

//#define SINGLE_FAST5_ONLY 1 //comment to enable multi fast5 support

// The following was adapted from Nanopolish : nanopolish_fast5_io.cpp

#ifndef SINGLE_FAST5_ONLY
#include <vector>
//
std::string fast5_get_raw_read_internal_name(fast5_file_t fh);
// Get a vector of read groups for a multi-fast5 file (eg [read_00041f-..., read_1243fe-....])
std::vector<std::string> fast5_get_multi_read_groups(fast5_file_t fh);
//
std::string fast5_get_string_attribute(fast5_file_t fh, const std::string& group_name, const std::string& attribute_name);
std::string fast5_get_raw_read_group(fast5_file_t fh, const std::string& read_id);
uint8_t fast5_is_vbz_compressed(fast5_file_t fh, const std::string& read_id);

#endif



#define LEGACY_FAST5_RAW_ROOT "/Raw/Reads/"


// from nanopolish_fast5_io.cpp
static inline fast5_file_t fast5_open(char* filename) {
    fast5_file_t fh;
    fh.hdf5_file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

#ifndef SINGLE_FAST5_ONLY
    // read and parse the file version to determine if this is a multi-fast5 structured file
    std::string version_str = fast5_get_string_attribute(fh, "/", "file_version");
    if(version_str != "") {
        int major;
        int minor;
        int ret = sscanf(version_str.c_str(), "%d.%d", &major, &minor);
        if(ret != 2) {
            fprintf(stderr, "Could not parse version string %s\n", version_str.c_str());
            exit(EXIT_FAILURE);
        }

        fh.is_multi_fast5 = major >= 1;
    } else {
        fh.is_multi_fast5 = false;
    }
#else
    fh.is_multi_fast5 = false;
#endif


    return fh;
}

//from nanopolish_fast5_io.cpp
static inline void fast5_close(fast5_file_t fh) {
    H5Fclose(fh.hdf5_file);
}

//from nanopolish_fast5_io.cpp
static inline float fast5_read_float_attribute(hid_t group,
                                               const char* attribute) {
    float val = NAN;

    hid_t attr = H5Aopen(group, attribute, H5P_DEFAULT);
    if (attr < 0) {
        WARNING("Failed to open attribute '%s'.", attribute);
        return NAN;
    }

    herr_t ret = H5Aread(attr, H5T_NATIVE_FLOAT, &val);
    if (ret < 0) {
        WARNING("Failed to read attribute '%s'.", attribute);
        return NAN;
    }

    H5Aclose(attr);

    return val;
}



static inline int32_t fast5_read_single_fast5(fast5_file_t fh, signal_t* sig) {

    hid_t hdf5_file = fh.hdf5_file;
    sig->rawptr = NULL;
    hid_t space;
    hsize_t nsample;
    herr_t status;

    // retrieve the size of the read name
    ssize_t size = H5Lget_name_by_idx(hdf5_file, LEGACY_FAST5_RAW_ROOT, H5_INDEX_NAME,
                                      H5_ITER_INC, 0, NULL, 0, H5P_DEFAULT);

    if (size < 0) {
        WARNING("Failed to retrieve the size of the read name.%s", "");
        return -1;
    }

    // copy the read name out of the fast5
    char* read_name = (char*)calloc(1 + size, sizeof(char));
    ssize_t ret =
        H5Lget_name_by_idx(hdf5_file, LEGACY_FAST5_RAW_ROOT, H5_INDEX_NAME, H5_ITER_INC, 0,
                           read_name, 1 + size, H5P_DEFAULT);
    if (ret < 0) {
        WARNING("Failed to retrieve the read name.%s", "");
        return -1;
    }

    // not the most efficient and safest, but being a bit lazy for the moment
    char* signal_path = (char*)malloc(
        (size + 1 + strlen(LEGACY_FAST5_RAW_ROOT) + strlen("/Signal")) * sizeof(char));
    MALLOC_CHK(signal_path);

    sprintf(signal_path, "%s%s%s", LEGACY_FAST5_RAW_ROOT, read_name, "/Signal");

#ifdef DEBUG_SIGNAL_PATH
    printf("Signal path : %s\n", signal_path);
#endif
    free(read_name);

    hid_t dset = H5Dopen(hdf5_file, signal_path, H5P_DEFAULT);
    if (dset < 0) {
        WARNING("Failed to open dataset '%s' to read raw signal.", signal_path);
        free(signal_path);
        return -1;
        // goto cleanup2;
    }

    space = H5Dget_space(dset);
    if (space < 0) {
        WARNING("Failed to create copy of dataspace for raw signal %s.",
                signal_path);
        H5Dclose(dset);
        free(signal_path);
        return -1;
        // goto cleanup3;
    }

    int32_t ret1 = H5Sget_simple_extent_dims(space, &nsample, NULL);
    if (ret1 < 0) {
        WARNING("Failed to get the dataspace dimension for raw signal %s.",
                signal_path);
        H5Sclose(space);
        H5Dclose(dset);
        free(signal_path);
        return -1;
    }

    sig->nsample = nsample;
    sig->rawptr = (float*)calloc(nsample, sizeof(float));
    status = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     sig->rawptr);

    if (status < 0) {
        free(sig->rawptr);
        WARNING("Failed to read raw data from dataset %s.", signal_path);
        H5Sclose(space);
        H5Dclose(dset);
        free(signal_path);
        return -1;
    }

    H5Sclose(space);
    H5Dclose(dset);

    // get channel parameters
    const char* scaling_path = "/UniqueGlobalKey/channel_id";

    hid_t scaling_group = H5Gopen(hdf5_file, scaling_path, H5P_DEFAULT);
    if (scaling_group < 0) {
        WARNING("Failed to open group %s.", scaling_path);
        free(signal_path);
        return -1;
    }

    sig->digitisation =
        fast5_read_float_attribute(scaling_group, "digitisation");
    sig->offset = fast5_read_float_attribute(scaling_group, "offset");
    sig->range = fast5_read_float_attribute(scaling_group, "range");
    sig->sample_rate =
        fast5_read_float_attribute(scaling_group, "sampling_rate");

    if (sig->digitisation == NAN || sig->offset == NAN || sig->range == NAN ||
        sig->sample_rate == NAN) {
        WARNING("Read a NAN value for scaling parameters for '%s'.",
                signal_path);
        H5Gclose(scaling_group);
        free(signal_path);
        return -1;
    }

    H5Gclose(scaling_group);
    free(signal_path);

    return 0;
}


#ifndef SINGLE_FAST5_ONLY

static inline int32_t fast5_read_multi_fast5(fast5_file_t fh, signal_t* sig, std::string read_id) {

    hid_t hdf5_file = fh.hdf5_file;
    sig->rawptr = NULL;
    hid_t space;
    hsize_t nsample;
    herr_t status;

    // mostly from scrappie
    std::string raw_read_group = fast5_get_raw_read_group(fh, read_id);

    // Create data set name
    std::string signal_path_str = raw_read_group + "/Signal";
    const char *signal_path = signal_path_str.c_str();

    hid_t dset = H5Dopen(hdf5_file, signal_path, H5P_DEFAULT);
    if (dset < 0) {
        WARNING("Failed to open dataset '%s' to read raw signal.", signal_path);
        return -1;
        // goto cleanup2;
    }

    space = H5Dget_space(dset);
    if (space < 0) {
        WARNING("Failed to create copy of dataspace for raw signal %s.",
                signal_path);
        H5Dclose(dset);
        return -1;
        // goto cleanup3;
    }

    int32_t ret1 = H5Sget_simple_extent_dims(space, &nsample, NULL);
    if (ret1 < 0) {
        WARNING("Failed to get the dataspace dimension for raw signal %s.",
                signal_path);
        H5Sclose(space);
        H5Dclose(dset);
        return -1;
    }

    sig->nsample = nsample;
    sig->rawptr = (float*)calloc(nsample, sizeof(float));
    status = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     sig->rawptr);

    if (status < 0) {
        free(sig->rawptr);
        if(fast5_is_vbz_compressed(fh, read_id) == 1) {
            ERROR("%s","The fast5 file is compressed with VBZ but the required plugin is not loaded. See https://f5c.page.link/troubleshoot for instructions.\n");
        }
        WARNING("Failed to read raw data from dataset %s.", signal_path);
        H5Sclose(space);
        H5Dclose(dset);
        return -1;
    }

    H5Sclose(space);
    H5Dclose(dset);

    // get channel parameters
    std::string scaling_path_str = fh.is_multi_fast5 ? "/read_" + read_id + "/channel_id"
                                                 :  "/UniqueGlobalKey/channel_id";
    const char *scaling_path = scaling_path_str.c_str();

    hid_t scaling_group = H5Gopen(hdf5_file, scaling_path, H5P_DEFAULT);
    if (scaling_group < 0) {
        WARNING("Failed to open group %s.", scaling_path);
        return -1;
    }

    sig->digitisation =
        fast5_read_float_attribute(scaling_group, "digitisation");
    sig->offset = fast5_read_float_attribute(scaling_group, "offset");
    sig->range = fast5_read_float_attribute(scaling_group, "range");
    sig->sample_rate =
        fast5_read_float_attribute(scaling_group, "sampling_rate");

    if (sig->digitisation == NAN || sig->offset == NAN || sig->range == NAN ||
        sig->sample_rate == NAN) {
        WARNING("Read a NAN value for scaling parameters for '%s'.",
                signal_path);
        H5Gclose(scaling_group);
        return -1;
    }

    H5Gclose(scaling_group);

    return 0;
}
#endif

static inline int32_t fast5_read(fast5_file_t fh, signal_t* sig, std::string read_id) {
#ifndef SINGLE_FAST5_ONLY
    if(fh.is_multi_fast5){
        return fast5_read_multi_fast5(fh, sig, read_id);
    }
    else{
        return fast5_read_single_fast5(fh, sig);
    }

#else
    return fast5_read_single_fast5(fh, sig);
#endif


}

//from nanopolish_fast5_io.cpp //used by the indexer
//todo convert to C
static inline std::string fast5_get_fixed_string_attribute(fast5_file_t fh, const std::string& group_name, const std::string& attribute_name)
{
    hid_t hdf5_file = fh.hdf5_file;
    size_t storage_size;
    char* buffer;
    hid_t group, attribute, attribute_type;
    int ret;
    std::string out;

    // according to http://hdf-forum.184993.n3.nabble.com/check-if-dataset-exists-td194725.html
    // we should use H5Lexists to check for the existence of a group/dataset using an arbitrary path
    ret = H5Lexists(hdf5_file, group_name.c_str(), H5P_DEFAULT);
    if(ret <= 0) {
        return "";
    }

    // Open the group /Raw/Reads/Read_nnn
    group = H5Gopen(hdf5_file, group_name.c_str(), H5P_DEFAULT);
    if(group < 0) {
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "could not open group %s\n", group_name.c_str());
#endif
        goto close_group;
    }

    // Ensure attribute exists
    ret = H5Aexists(group, attribute_name.c_str());
    if(ret <= 0) {
        goto close_group;
    }

    // Open the attribute
    attribute = H5Aopen(group, attribute_name.c_str(), H5P_DEFAULT);
    if(attribute < 0) {
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "could not open attribute: %s\n", attribute_name.c_str());
#endif
        goto close_attr;
    }

    // Get data type and check it is a fixed-length string
    attribute_type = H5Aget_type(attribute);
    if(H5Tis_variable_str(attribute_type)) {
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "variable length string detected -- ignoring attribute\n");
#endif
        goto close_type;
    }

    // Get the storage size and allocate
    storage_size = H5Aget_storage_size(attribute);
    buffer = (char*)calloc(storage_size + 1, sizeof(char));

    // finally read the attribute
    ret = H5Aread(attribute, attribute_type, buffer);
    if(ret >= 0) {
        out = buffer;
    }

    // clean up
    free(buffer);
close_type:
    H5Tclose(attribute_type);
close_attr:
    H5Aclose(attribute);
close_group:
    H5Gclose(group);

    return out;
}

static inline std::string fast5_get_read_id_single_fast5(fast5_file_t fh)
{

    // this function is not supported for multi-fast5 files
    assert(!fh.is_multi_fast5);

    hid_t hdf5_file = fh.hdf5_file;
    std::string out = "";

    // Get the path to the raw read group
    // retrieve the size of the read name
    ssize_t size = H5Lget_name_by_idx(hdf5_file, LEGACY_FAST5_RAW_ROOT, H5_INDEX_NAME,
                                      H5_ITER_INC, 0, NULL, 0, H5P_DEFAULT);

    if (size < 0) {
        WARNING("Failed to retrieve the size of the read name.%s", "");
        return out;
    }

    // copy the read name out of the fast5
    char* read_name = (char*)calloc(1 + size, sizeof(char));
    ssize_t ret =
        H5Lget_name_by_idx(hdf5_file, LEGACY_FAST5_RAW_ROOT, H5_INDEX_NAME, H5_ITER_INC, 0,
                           read_name, 1 + size, H5P_DEFAULT);
    if (ret < 0) {
        WARNING("Failed to retrieve the read name.%s", "");
        return out;
    }

    std::string raw_read_group= std::string(LEGACY_FAST5_RAW_ROOT) + read_name;
    free(read_name);


    return fast5_get_fixed_string_attribute(fh, raw_read_group, "read_id");
}

#else

static inline fast5_file_t fast5_open(char* filename) {
    ERROR("%s", "f5c has been compiled with no FAST5/HDF5 support. s2f unavailable. Recompile with FAST5/HDF5 support.");
    exit(EXIT_FAILURE);
}

static inline int32_t fast5_read(fast5_file_t fh, signal_t* sig, std::string read_id){
    ERROR("%s", "f5c has been compiled with no FAST5/HDF5 support. s2f unavailable. Recompile with FAST5/HDF5 support.");
    exit(EXIT_FAILURE);
}

static inline void fast5_close(fast5_file_t fh) {
    ERROR("%s", "f5c has been compiled with no FAST5/HDF5 support. s2f unavailable. Recompile with FAST5/HDF5 support.");
    exit(EXIT_FAILURE);
}

static inline std::string fast5_get_read_id_single_fast5(fast5_file_t fh){
    ERROR("%s", "f5c has been compiled with no FAST5/HDF5 support. s2f unavailable. Recompile with FAST5/HDF5 support.");
    exit(EXIT_FAILURE);
}

#include <vector>
static inline std::vector<std::string> fast5_get_multi_read_groups(fast5_file_t fh){
    ERROR("%s", "f5c has been compiled with no FAST5/HDF5 support. s2f unavailable. Recompile with FAST5/HDF5 support.");
    exit(EXIT_FAILURE);
}

#endif

#endif
