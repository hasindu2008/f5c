/* @file index.c
**
** implementation of f5c index
** Code was copied/adapted from Nanopolish index module originally authored by Jared Simpson
** Code was extended by Hasindu Gamaarachchi to perform indexing in parallel using multiple processes
** @@
******************************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <getopt.h>
#include <assert.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include <slow5/slow5.h>

#include "nanopolish_read_db.h"
#include "error.h"
#include <vector>
#include <dirent.h>
#include "f5c.h"
#include "fast5lite.h"
#include "f5cmisc.h"

// from nanopolish
// ref: http://stackoverflow.com/a/612176/717706
// return true if the given name is a directory
bool is_directory(const std::string& file_name){
    auto dir = opendir(file_name.c_str());
    if(not dir) {
        return false;
    }
    closedir(dir);
    return true;
}

// from nanopolish
// Split a string into parts based on the delimiter (nanopolish_common.cpp)
std::vector<std::string> split(std::string in, char delimiter);

// from nanopolish
std::vector< std::string > list_directory(const std::string& file_name)
{
    std::vector< std::string > res;
    DIR* dir;
    struct dirent *ent;

    dir = opendir(file_name.c_str());
    if(not dir) {
        return res;
    }
    while((ent = readdir(dir)) != nullptr) {
        res.push_back(ent->d_name);
    }
    closedir(dir);
    return res;
}


static const char *INDEX_USAGE_MESSAGE =
"Usage: f5c index [OPTIONS] -d fast5_directory reads.fastq\n"
"       f5c index [OPTIONS] --slow5 signals.blow5 reads.fastq\n\n"
"Build an index for accessing the base sequence (fastq/fasta) and raw signal (fast5/slow5) for a given read ID. f5c index is an extended and optimised version of nanopolish index by Jared Simpson\n"
"\n"
"  -h                display this help and exit\n"
"  -d STR            path to the directory containing fast5 files. This option can be given multiple times.\n"
"  -s STR            the sequencing summary file\n"
"  -f STR            file containing the paths to the sequencing summary files (one per line)\n"
"  -t INT            number of threads used for bgzf compression (makes indexing faster)\n"
"  --iop INT         number of I/O processes to read fast5 files (makes fast5 indexing faster)\n"
"  --slow5 FILE      slow5 file containing raw signals\n"
"  --skip-slow5-idx  do not build the .idx for the slow5 file (useful when a slow5 index is already available)\n"
"  --verbose INT     verbosity level\n"
"  --version         print version\n"
"\nSee the manual page for details (`man ./docs/f5c.1' or https://f5c.page.link/man)."
"\n\n"
;

namespace opt
{
    static unsigned int verbose = 1;
    static std::vector<std::string> raw_file_directories;
    static std::string reads_file;
    static std::vector<std::string> sequencing_summary_files;
    static std::string sequencing_summary_fofn;
    int iop = 1;
    int threads = 4;
    char *slow5file = NULL;
    bool skip_slow5_idx = false;
}
//static std::ostream* os_p;

// from nanopolish
void index_file_from_map(ReadDB& read_db, const std::string& fn, const std::multimap<std::string, std::string>& fast5_to_read_name_map)
{
    //PROFILE_FUNC("index_file_from_map")

    // Check if this fast5 file is in the map
    size_t last_dir_pos = fn.find_last_of("/");
    std::string fast5_basename = last_dir_pos == std::string::npos ? fn : fn.substr(last_dir_pos + 1);

    auto range = fast5_to_read_name_map.equal_range(fast5_basename);
    for(auto iter = range.first; iter != range.second; ++iter) {
        if(read_db.has_read(iter->second)) {
            read_db.add_signal_path(iter->second.c_str(), fn);
        }
    }
}

// adapted from nanopolish_index
void index_file_from_fast5(ReadDB& read_db, const std::string& fn)
{
    //PROFILE_FUNC("index_file_from_fast5")

    char* fast5_path =(char*)malloc(fn.size()+1);
    strcpy(fast5_path, fn.c_str());

    fast5_file_t f5_file =  fast5_open(fast5_path);
    hid_t hdf5_file = f5_file.hdf5_file;
    if(hdf5_file < 0) {
        fprintf(stderr, "could not open fast5 file: %s\n", fast5_path);
    }

    if(f5_file.is_multi_fast5) {
        std::vector<std::string> read_groups = fast5_get_multi_read_groups(f5_file);
        std::string prefix = "read_";
        for(size_t group_idx = 0; group_idx < read_groups.size(); ++group_idx) {
            std::string group_name = read_groups[group_idx];
            if(group_name.find(prefix) == 0) {
                std::string read_id = group_name.substr(prefix.size());
                read_db.add_signal_path(read_id, fn);
            }
        }
    } else {
        std::string read_id = fast5_get_read_id_single_fast5(f5_file);
        if(read_id != "") {
            read_db.add_signal_path(read_id, fn);
        }
    }
    free(fast5_path);
    fast5_close(f5_file);
}

// from nanopolish
void index_path(ReadDB& read_db, const std::string& path, const std::multimap<std::string, std::string>& fast5_to_read_name_map)
{
    fprintf(stderr, "[readdb] indexing %s\n", path.c_str());
    if (is_directory(path)) {
        auto dir_list = list_directory(path);
        for (const auto& fn : dir_list) {
            if(fn == "." or fn == "..") {
                continue;
            }

            std::string full_fn = path + "/" + fn;
            bool is_fast5 = full_fn.find(".fast5") != std::string::npos;
            bool in_map = fast5_to_read_name_map.find(fn) != fast5_to_read_name_map.end();

            // JTS 04/19: is_directory is painfully slow so we first check if the file is in the name map
            // if it is, it is definitely not a directory so we can skip the system call
            if(!in_map && is_directory(full_fn)) {
                // recurse
                index_path(read_db, full_fn, fast5_to_read_name_map);
            } else if (is_fast5) {
                if(in_map) {
                    index_file_from_map(read_db, full_fn, fast5_to_read_name_map);
                } else {
                    index_file_from_fast5(read_db, full_fn);
                }
            }
        }
    }
}

// from nanopolish
// read sequencing summary files from the fofn and add them to the list
void process_summary_fofn()
{
    if(opt::sequencing_summary_fofn.empty()) {
        return;
    }

    // open
    std::ifstream in_file(opt::sequencing_summary_fofn.c_str());
    if(in_file.bad()) {
        fprintf(stderr, "error: could not file %s\n", opt::sequencing_summary_fofn.c_str());
        exit(EXIT_FAILURE);
    }

    // read
    std::string filename;
    while(getline(in_file, filename)) {
        opt::sequencing_summary_files.push_back(filename);
    }
}

// from nanopolish
void exit_bad_header(const std::string& str, const std::string& filename)
{
    fprintf(stderr, "Could not find %s column in the header of %s\n", str.c_str(), filename.c_str());
    exit(EXIT_FAILURE);
}

// from nanopolish
void parse_sequencing_summary(const std::string& filename, std::multimap<std::string, std::string>& out_map)
{
    // open
    std::ifstream in_file(filename.c_str());
    if(!in_file.good()) {
        fprintf(stderr, "error: could not read file %s\n", filename.c_str());
        exit(EXIT_FAILURE);
    }

    // read header to get the column index of the read and file name
    std::string header;
    getline(in_file, header);
    std::vector<std::string> fields = split(header, '\t');

    const std::string READ_NAME_STR = "read_id";
    size_t filename_idx = SIZE_MAX;
    size_t read_name_idx = SIZE_MAX;

    for(size_t i = 0; i < fields.size(); ++i) {
        if(fields[i] == READ_NAME_STR) {
            read_name_idx = i;
        }

        // 19/11/05: support live basecalling summary files
        if(fields[i] == "filename" || fields[i] == "filename_fast5") {
            filename_idx = i;
        }
    }

    if(filename_idx == SIZE_MAX ) {
        exit_bad_header("fast5 filename", filename);
    }

    if(read_name_idx == SIZE_MAX ) {
        exit_bad_header(READ_NAME_STR, filename);
    }

    // read records and add to map
    std::string line;
    while(getline(in_file, line)) {
        fields = split(line, '\t');
        std::string fast5_filename = fields[filename_idx];
        std::string read_name = fields[read_name_idx];
        out_map.insert(std::make_pair(fast5_filename, read_name));
    }
}

static const char* shortopts = "Vhd:f:s:v:t:";

static const struct option longopts[] = {
    { "help",                      no_argument,       NULL, 'h' },     //0
    { "version",                   no_argument,       NULL, 'V' },  //1
    { "verbose",                   no_argument,       NULL, 'v' },          //2
    { "directory",                 required_argument, NULL, 'd' },          //3
    { "sequencing-summary-file",   required_argument, NULL, 's' },          //4
    { "summary-fofn",              required_argument, NULL, 'f' },          //5
    { "iop",                       required_argument, NULL, 0},             //6
    { "threads",                   required_argument, NULL, 't'},           //7
    { "slow5",                     required_argument, NULL, 0},           //8
    { "skip-slow5-idx",            no_argument,       NULL, 0},           //9
    { NULL, 0, NULL, 0 }
};

// adapted from nanopolish
void parse_index_options(int argc, char** argv)
{
    bool die = false;
    int longindex = 0;
    //std::vector< std::string> log_level;
    for (int c; (c = getopt_long(argc, argv, shortopts, longopts, &longindex)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'h':
                std::cout << INDEX_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case 'V':
                fprintf(stdout,"F5C %s\n",F5C_VERSION);
                exit(EXIT_SUCCESS);
            case 'v': arg >> opt::verbose; break;
            case 's': opt::sequencing_summary_files.push_back(arg.str()); break;
            case 'd': opt::raw_file_directories.push_back(arg.str()); break;
            case 'f': arg >> opt::sequencing_summary_fofn; break;
            case 't':
                arg >> opt::threads;
                if (opt::threads < 1) {
                    ERROR("Number of threads should be larger than 0. You entered %d", opt::threads);
                    exit(EXIT_FAILURE);
                }
                break;
            case  0 :
                if (longindex == 6) {
                    arg >> opt::iop;
                    if (opt::iop < 1) {
                        ERROR("Number of I/O processes should be larger than 0. You entered %d", opt::iop);
                        exit(EXIT_FAILURE);
                    }
                }
                if (longindex == 8) {
                    opt::slow5file = optarg;
                }
                if (longindex == 9) {
                    opt::skip_slow5_idx = true;
                }
                break;
        }
    }

    if (argc - optind < 1) {
        std::cerr << "[f5c index]  not enough arguments\n";
        die = true;
    }

    if (argc - optind > 1) {
        std::cerr << "[f5c index] too many arguments\n";
        die = true;
    }

    if (die)
    {
        std::cout << "\n" << INDEX_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    if(opt::slow5file == NULL){
        if(opt::iop > 1){
            if(!opt::sequencing_summary_fofn.empty()){
                WARNING("%s","--iop is incompatible with sequencing summary files. Option --summary-fofn will be ignored");
            }
            if(!opt::sequencing_summary_files.empty()){
                WARNING("%s","--iop is incompatible with sequencing summary files. Option --sequencing-summary-file will be ignored");
            }
        }
        if(opt::skip_slow5_idx){
            WARNING("%s","--skip-slow5-idx is only useful with --slow5. Option ignored");
        }
        // else {
        //     if(opt::sequencing_summary_fofn.empty() && opt::sequencing_summary_files.empty()){
        //         INFO("%s","Consider using --iop option for fast parallel indexing");
        //     }
        // }
        INFO("%s","Consider using --slow5 option for fast indexing, methylation calling and eventalignment. See f5c section under https://hasindu2008.github.io/slow5tools/workflows.html for an example.");
    } else{
        if(opt::iop > 1){
            WARNING("%s","Option --iop has no effect in slow5 mode");
        }
        if(!opt::sequencing_summary_fofn.empty()){
            WARNING("%s","Option --summary-fofn will be ignored in slow5 mode");
        }
        if(!opt::sequencing_summary_files.empty()){
            WARNING("%s","Option --sequencing-summary-file will be ignored in slow5 mode");
        }
        if(!opt::raw_file_directories.empty()){
            WARNING("%s","Option -d will be ignored in slow5 mode");
        }
    }

    if(opt::threads == 1){
        INFO("%s","Consider using -t option for fast multi-threaded bgzf");
    }

    opt::reads_file = argv[optind++];
}

// from nanopolish
void clean_fast5_map(std::multimap<std::string, std::string>& mmap)
{
    std::map<std::string, int> fast5_read_count;
    for(const auto& iter : mmap) {
        fast5_read_count[iter.first]++;
    }

    // calculate the mode of the read counts per fast5
    std::map<size_t, size_t> read_count_distribution;
    for(const auto& iter : fast5_read_count) {
        read_count_distribution[iter.second]++;
    }

    size_t mode = 0;
    size_t mode_count = 0;
    for(const auto& iter : read_count_distribution) {
        if(iter.second > mode_count) {
            mode = iter.first;
            mode_count = iter.second;
        }
    }

    // this is the default number of reads per fast5 and we expect the summary file to support that
    // this code is to detect whether the user selected a different number of reads (ARTIC recommends 1000)
    // so the summary file still works
    int EXPECTED_ENTRIES_PER_FAST5 = 4000;
    if(mode_count > 20) {
        EXPECTED_ENTRIES_PER_FAST5 = mode;
    }

    std::vector<std::string> invalid_entries;
    for(const auto& iter : fast5_read_count) {
        if(iter.second > EXPECTED_ENTRIES_PER_FAST5) {
            //fprintf(stderr, "warning: %s has %d entries in the summary and will be indexed the slow way\n", iter.first.c_str(), iter.second);
            invalid_entries.push_back(iter.first);
        }
    }

    if(invalid_entries.size() > 0) {
        fprintf(stderr, "warning: detected invalid summary file entries for %zu of %zu fast5 files\n", invalid_entries.size(), fast5_read_count.size());
        fprintf(stderr, "These files will be indexed without using the summary file, which is slow.\n");
    }

    for(const auto& fast5_name : invalid_entries) {
        mmap.erase(fast5_name);
    }
}



// given a directory path, recursively find all fast5 files
void find_all_fast5(const std::string& path, std::vector<std::string>& fast5_files)
{
    STDERR("Looking for fast5 in %s", path.c_str());
    if (is_directory(path)) {
        auto dir_list = list_directory(path);
        for (const auto& fn : dir_list) {
            if(fn == "." or fn == "..") {
                continue;
            }

            std::string full_fn = path + "/" + fn;
            bool is_fast5 = full_fn.find(".fast5") != std::string::npos;
            // JTS 04/19: is_directory is painfully slow
            if(is_directory(full_fn)) {
                // recurse
                find_all_fast5(full_fn,fast5_files);
            } else if (is_fast5) {
                fast5_files.push_back(full_fn);
                //add to the list
            }
        }
    }
}

// args for processes
typedef struct {
    int32_t starti;
    int32_t endi;
    int32_t proc_index;
    std::string tmp_file;
}proc_arg_t;


// find the readids in a given fast5 file and dump to the tmp_file
// todo : repeated content with index_file_from_fast5, can be modularised
void get_readids_from_fast5(const std::string& fn, FILE *tmp_file)
{
    //PROFILE_FUNC("index_file_from_fast5")

    char* fast5_path =(char*)malloc(fn.size()+1);
    strcpy(fast5_path, fn.c_str());

    fast5_file_t f5_file =  fast5_open(fast5_path);
    hid_t hdf5_file = f5_file.hdf5_file;
    if(hdf5_file < 0) {
        fprintf(stderr, "could not open fast5 file: %s\n", fast5_path);
    }

    if(f5_file.is_multi_fast5) {
        std::vector<std::string> read_groups = fast5_get_multi_read_groups(f5_file);
        std::string prefix = "read_";
        for(size_t group_idx = 0; group_idx < read_groups.size(); ++group_idx) {
            std::string group_name = read_groups[group_idx];
            if(group_name.find(prefix) == 0) {
                std::string read_id = group_name.substr(prefix.size());
                fprintf(tmp_file,"%s\t%s\n",read_id.c_str(), fast5_path);
            }
        }
    } else {
        std::string read_id = fast5_get_read_id_single_fast5(f5_file);
        if(read_id != "") {
            fprintf(tmp_file,"%s\t%s\n",read_id.c_str(), fast5_path);
        }
    }
    free(fast5_path);
    fast5_close(f5_file);
}


// what a child process should do, i.e. open a tmp file, go through the fast5 files
void f5c_index_child_worker(proc_arg_t args, const std::vector<std::string>& fast5_files){

    FILE *tmp_file = fopen(args.tmp_file.c_str(),"w");
    F_CHK(tmp_file,args.tmp_file.c_str());

    int i=0;
    for (i = args.starti; i < args.endi; i++) {
        get_readids_from_fast5(fast5_files[i], tmp_file);
    }

    fclose(tmp_file);

}

// find all fast5 files in the specified directories, divide the work and spawn multiple I/O processes
void f5c_index_iop(int iop, std::string reads_file, std::vector<std::string> raw_file_directories, int verbosity){

    double realtime0 = realtime();
    std::vector<std::string> fast5_files;
    for(const auto& dir_name : opt::raw_file_directories) {
        find_all_fast5(dir_name, fast5_files);
    }
    int64_t num_fast5_files = fast5_files.size();
    fprintf(stderr, "[%s] %ld fast5 files found - took %.3fs\n", __func__, fast5_files.size(), realtime() - realtime0);


    realtime0 = realtime();
    //create processes
    pid_t pids[iop];
    proc_arg_t proc_args[iop];
    int32_t t;
    int32_t i = 0;
    int32_t step = (num_fast5_files + iop - 1) / iop;
    //todo : check for higher num of procs than the data
    //current works but many procs are created despite

    //set the data structures
    for (t = 0; t < iop; t++) {
        proc_args[t].starti = i;
        i += step;
        if (i > num_fast5_files) {
            proc_args[t].endi = num_fast5_files;
        } else {
            proc_args[t].endi = i;
        }
        proc_args[t].proc_index = t;

        char tmp[11];
        sprintf(tmp, "%d",t);
        proc_args[t].tmp_file = reads_file + ".readdb.tmp" + tmp;
    }

    //create processes
    STDERR("Spawning %d I/O processes to circumvent HDF hell",iop);
    for(t = 0; t < iop; t++){
        pids[t] = fork();

        if(pids[t]==-1){
            ERROR("%s","Fork failed");
            perror("");
            exit(EXIT_FAILURE);
        }
        if(pids[t]==0){ //child
            f5c_index_child_worker(proc_args[t],fast5_files);
            exit(EXIT_SUCCESS);
        }
        if(pids[t]>0){ //parent
            continue;
        }
    }

    //wait for processes
    int status,w;
    for (t = 0; t < iop; t++) {
        if(verbosity>1){
            STDERR("parent : Waiting for child with pid %d",pids[t]);
        }

            w = waitpid(pids[t], &status, 0);
            if (w == -1) {
                ERROR("%s","waitpid failed");
                perror("");
                exit(EXIT_FAILURE);
            }
            else if (WIFEXITED(status)){
                if(verbosity>1){
                    STDERR("child process %d exited, status=%d", pids[t], WEXITSTATUS(status));
                }
                if(WEXITSTATUS(status)!=0){
                    ERROR("child process %d exited with status=%d",pids[t], WEXITSTATUS(status));
                    exit(EXIT_FAILURE);
                }
            }
            else {
                if (WIFSIGNALED(status)) {
                    ERROR("child process %d killed by signal %d", pids[t], WTERMSIG(status));
                } else if (WIFSTOPPED(status)) {
                    ERROR("child process %d stopped by signal %d", pids[t], WSTOPSIG(status));
                } else {
                    ERROR("child process %d did not exit propoerly: status %d", pids[t], status);
                }
                exit(EXIT_FAILURE);
            }


    }
    fprintf(stderr, "[%s] Parallel indexing done - took %.3fs\n", __func__,  realtime() - realtime0);

}

// merge temporary dumps produced by child processes
void f5c_index_merge(ReadDB& read_db, int iop, std::string reads_file){

    double realtime0 = realtime();
    int i;
    for(i=0;i<iop;i++){
        char tmp[11];
        sprintf(tmp, "%d",i);
        std::string filename = reads_file + ".readdb.tmp" + tmp;

        std::ifstream in_file(filename.c_str());
        if(!in_file.good()) {
            ERROR("could not read the temporary file %s\n", filename.c_str());
            exit(EXIT_FAILURE);
        }

        // read the database
        std::string line;
        while(getline(in_file, line)) {
            std::vector<std::string> fields = split(line, '\t');

            std::string name = "";
            std::string path = "";
            if(fields.size() == 2) {
                name = fields[0];
                path = fields[1];
                read_db.add_signal_path(name, path);
            }
        }
        in_file.close();
        int ret = remove(filename.c_str());
        if(ret != 0){
            WARNING("Removing temporary file %s failed",filename.c_str());
            perror("");
        }

    }
    fprintf(stderr, "[%s] Indexing merging done - took %.3fs.\n", __func__, realtime() - realtime0);

}

// adapted from nanopolish
int index_main(int argc, char** argv)
{
    parse_index_options(argc, argv);

    if(opt::slow5file == NULL){

        std::multimap<std::string, std::string> fast5_to_read_name_map;

        if(opt::iop == 1){

            // Read a map from fast5 files to read name from the sequencing summary files (if any)
            process_summary_fofn();
            for(const auto& ss_filename : opt::sequencing_summary_files) {
                if(opt::verbose > 2) {
                    fprintf(stderr, "summary: %s\n", ss_filename.c_str());
                }
                parse_sequencing_summary(ss_filename, fast5_to_read_name_map);
            }

            // Detect non-unique fast5 file names in the summary file
            // This occurs when a file with the same name (abc_0.fast5) appears in both fast5_pass
            // and fast5_fail. This will be fixed by ONT but in the meantime we check for
            // fast5 files that have a non-standard number of reads (>4000) and remove them
            // from the map. These fast5s will be indexed the slow way.
            clean_fast5_map(fast5_to_read_name_map);

        }

        if(opt::iop>1){ //forking is better before doing big memory allocations
            f5c_index_iop(opt::iop, opt::reads_file, opt::raw_file_directories, opt::verbose);
        }

        // import read names, and possibly fast5 paths, from the fasta/fastq file
        double realtime0 = realtime();
        ReadDB read_db;
        read_db.build(opt::reads_file, opt::threads);
        bool all_reads_have_paths = read_db.check_signal_paths();
        if(opt::verbose > 1) fprintf(stderr, "[%s] Readdb built - took %.3fs\n", __func__, realtime() - realtime0);

        // if the input fastq did not contain a complete set of paths
        // use the fofn/directory provided to augment the index
        if(opt::iop==1){
            realtime0 = realtime();
            if(!all_reads_have_paths) {
                for(const auto& dir_name : opt::raw_file_directories) {
                    index_path(read_db, dir_name, fast5_to_read_name_map);
                }
            }
            if(opt::verbose > 1) fprintf(stderr, "[%s] Indexing done - took %.3fs\n", __func__, realtime() - realtime0);
        }
        else{
            f5c_index_merge(read_db, opt::iop, opt::reads_file);
        }

        realtime0 = realtime();
        size_t num_with_path = read_db.get_num_reads_with_path();
        size_t num_reads = read_db.get_num_reads();
        if(num_with_path == 0) {
            ERROR("%s","No fast5 files found");
            exit(EXIT_FAILURE);
        } else {
            read_db.print_stats();
            read_db.save();
            if (num_with_path < num_reads * 0.99 ) {
                WARNING("fast5 files could not be located for %ld reads",num_reads-num_with_path);
            }
        }
        if(opt::verbose > 1) fprintf(stderr, "[%s] Readdb saved - took %.3fs\n", __func__, realtime() - realtime0);
    }
    else{

        double realtime0 = realtime();
        slow5_file_t *sp = slow5_open(opt::slow5file,"r");
        if(sp==NULL){
            ERROR("Error in opening slowfile %s",opt::slow5file);
            exit(EXIT_FAILURE);
        }
        if(!opt::skip_slow5_idx){
            int ret=0;
            ret = slow5_idx_create(sp);
            if(ret<0){
                ERROR("Error in creating index for slow5 file %s\n",opt::slow5file);
                exit(EXIT_FAILURE);
            }
            if(opt::verbose > 0) fprintf(stderr, "[%s] Slow5 index built - took %.3fs\n", __func__, realtime() - realtime0);
        }
        slow5_close(sp);

        realtime0 = realtime();
        ReadDB read_db;
        read_db.build(opt::reads_file, opt::threads);
        read_db.clean();
        if(opt::verbose > 0) fprintf(stderr, "[%s] Fasta index built - took %.3fs\n", __func__, realtime() - realtime0);

    }

    return 0;
}
