//---------------------------------------------------------
// Copyright 2018 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_fast5_io -- lightweight functions
// to read specific data from fast5 files
//
#include <string.h>
#include <math.h>
#include <assert.h>
#include "f5c.h"
#include "fast5lite.h"

#ifndef DISABLE_HDF5

#ifndef SINGLE_FAST5_ONLY

#define H5Z_FILTER_VBZ 32020 //We need to find out what the numerical value for this is

uint8_t fast5_is_vbz_compressed(fast5_file_t fh, const std::string& read_id) {

    hid_t dset, dcpl;
    H5Z_filter_t filter_id = 0;
    char filter_name[80];
    size_t nelmts = 1; /* number of elements in cd_values */
    unsigned int values_out[1] = {99};
    unsigned int flags;

    // mostly from scrappie
    std::string raw_read_group = fast5_get_raw_read_group(fh, read_id);

    // Create data set name
    std::string signal_path = raw_read_group + "/Signal";

    dset = H5Dopen (fh.hdf5_file, signal_path.c_str(), H5P_DEFAULT);

    dcpl = H5Dget_create_plist (dset);

    filter_id = H5Pget_filter2 (dcpl, (unsigned) 0, &flags, &nelmts, values_out, sizeof(filter_name) - 1, filter_name, NULL);

    H5Pclose (dcpl);
    H5Dclose (dset);

    if(filter_id == H5Z_FILTER_VBZ){
        return 1;
    }
    else{
        return 0;
    }
}

//
bool fast5_is_open(fast5_file_t fh)
{
    return fh.hdf5_file >= 0;
}

//
std::vector<std::string> fast5_get_multi_read_groups(fast5_file_t fh)
{
    std::vector<std::string> out;
    ssize_t buffer_size = 0;
    char* buffer = NULL;

    // get the number of groups in the root group
    H5G_info_t group_info;
    int ret = H5Gget_info_by_name(fh.hdf5_file, "/", &group_info, H5P_DEFAULT);
    if(ret < 0) {
        fprintf(stderr, "error getting group info\n");
        exit(EXIT_FAILURE);
    }

    for(size_t group_idx = 0; group_idx < group_info.nlinks; ++group_idx) {

        // retrieve the size of this group name
        ssize_t size = H5Lget_name_by_idx(fh.hdf5_file, "/", H5_INDEX_NAME, H5_ITER_INC, group_idx, NULL, 0, H5P_DEFAULT);

        if(size < 0) {
            fprintf(stderr, "error getting group name size\n");
            exit(EXIT_FAILURE);
        }
        size += 1; // for null terminator

        if(size > buffer_size) {
            buffer = (char*)realloc(buffer, size);
            buffer_size = size;
        }

        // copy the group name
        H5Lget_name_by_idx(fh.hdf5_file, "/", H5_INDEX_NAME, H5_ITER_INC, group_idx, buffer, buffer_size, H5P_DEFAULT);
        buffer[size-1] = '\0';
        out.push_back(buffer);
    }

    free(buffer);
    buffer = NULL;
    buffer_size = 0;
    return out;
}


std::string fast5_get_experiment_type(fast5_file_t fh, const std::string& read_id)
{
    std::string group = fh.is_multi_fast5 ? "/read_" + read_id + "/context_tags"
                                          : "/UniqueGlobalKey/context_tags";
    return fast5_get_string_attribute(fh, group.c_str(), "experiment_type");
}

//
// Internal functions
//

//
std::string fast5_get_raw_root(fast5_file_t fh, const std::string& read_id)
{
    return fh.is_multi_fast5 ? "/read_" + read_id + "/Raw" : "/Raw/Reads";
}

//
std::string fast5_get_raw_read_group(fast5_file_t fh, const std::string& read_id)
{
    if(fh.is_multi_fast5) {
        return "/read_" + read_id + "/Raw";
    } else {
        std::string internal_read_name = fast5_get_raw_read_internal_name(fh);
        return internal_read_name != "" ? std::string(LEGACY_FAST5_RAW_ROOT) + "/" + internal_read_name : "";
    }
}

//
std::string fast5_get_raw_read_internal_name(fast5_file_t fh)
{
    // This code is From scrappie's fast5_interface

    // retrieve the size of the read name
    ssize_t size =
        H5Lget_name_by_idx(fh.hdf5_file, LEGACY_FAST5_RAW_ROOT, H5_INDEX_NAME, H5_ITER_INC, 0, NULL, 0, H5P_DEFAULT);

    if (size < 0) {
        return "";
    }

    // copy the read name out of the fast5
    char* name = (char*)calloc(1 + size, sizeof(char));
    H5Lget_name_by_idx(fh.hdf5_file, LEGACY_FAST5_RAW_ROOT, H5_INDEX_NAME, H5_ITER_INC, 0, name, 1 + size, H5P_DEFAULT);

    // cleanup
    std::string out(name);
    free(name);
    return out;
}

//
std::string fast5_get_string_attribute(fast5_file_t fh, const std::string& group_name, const std::string& attribute_name)
{
    hid_t group, attribute, attribute_type, native_type;
    std::string out;

    // according to http://hdf-forum.184993.n3.nabble.com/check-if-dataset-exists-td194725.html
    // we should use H5Lexists to check for the existence of a group/dataset using an arbitrary path
    // HDF5 1.8 returns 0 on the root group, so we explicitly check for it
    int ret = group_name == "/" ? 1 : H5Lexists(fh.hdf5_file, group_name.c_str(), H5P_DEFAULT);
    if(ret <= 0) {
        return "";
    }

    // Open the group containing the attribute we want
    group = H5Gopen(fh.hdf5_file, group_name.c_str(), H5P_DEFAULT);
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
    if(attribute_type < 0) {
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "failed to get attribute type %s\n", attribute_name.c_str());
#endif
        goto close_type;
    }

    if(H5Tget_class(attribute_type) != H5T_STRING) {
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "attribute %s is not a string\n", attribute_name.c_str());
#endif
        goto close_type;
    }

    native_type = H5Tget_native_type(attribute_type, H5T_DIR_ASCEND);
    if(native_type < 0) {
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "failed to get native type for %s\n", attribute_name.c_str());
#endif
        goto close_native_type;
    }

    if(H5Tis_variable_str(attribute_type) > 0) {
        // variable length string
        char* buffer;
        ret = H5Aread(attribute, native_type, &buffer);
        if(ret < 0) {
            fprintf(stderr, "error reading attribute %s\n", attribute_name.c_str());
            exit(EXIT_FAILURE);
        }
        out = buffer;
        free(buffer);
        buffer = NULL;

    } else {
        // fixed length string
        size_t storage_size;
        char* buffer;

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
    }

close_native_type:
    H5Tclose(native_type);
close_type:
    H5Tclose(attribute_type);
close_attr:
    H5Aclose(attribute);
close_group:
    H5Gclose(group);

    return out;
}

#endif

#endif
