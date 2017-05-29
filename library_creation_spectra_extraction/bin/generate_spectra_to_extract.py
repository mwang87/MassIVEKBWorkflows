#!/usr/bin/python


import sys
import getopt
import os
import json
import ming_psm_library
import ming_fileio_library
import ming_proteosafe_library
import library_creation
import pickle
from collections import defaultdict

def usage():
    output_usage = "<parallel json> <input params> <input folder of results>"
    output_usage += " <output folder>"
    print(output_usage)

def main():
    parallel_json = json.loads(open(sys.argv[1]).read())
    params_filename = sys.argv[2]
    input_folder_of_results = sys.argv[3]
    output_folder = sys.argv[4]

    my_node = parallel_json["node_partition"]
    total_node = parallel_json["total_paritions"]

    all_input_files = ming_fileio_library.list_files_in_dir(input_folder_of_results)
    all_input_files.sort()

    ###
    ### TODO We will have to read parameters and see if we need to eliminate some PSMs, with PSM FDR filter, KL Filter, ambiguity score filter, unique intensity filter
    ###

    params_obj = ming_proteosafe_library.parse_xml_file(open(params_filename))
    total_file_count = 0
    all_input_files = all_input_files[my_node::total_node]
    current_working_psm_set = ming_psm_library.PSMset("Ming")

    for input_file in all_input_files:
        #Assume these are variant files
        #We can treat this like a psm file and then combine all of the as a new variants file
        total_file_count += 1
        print(input_file, total_file_count, "of", len(all_input_files))
        input_pickle = open(input_file, 'rb')
        temp_psm_set = pickle.load(input_pickle)
        print("Loaded", len(temp_psm_set.psms))

        for psm in temp_psm_set.psms:
            precursor_string = "%s:%d" % (psm.annotation, psm.charge)
            score = psm.score

            #Determine minimum score cutoff
            current_score = psm.sorting_value()
            peptide_length = len(psm.get_stripped_sequence())

            current_working_psm_set.psms.append(psm)

    #Saving out psms
    output_filename = os.path.join(output_folder, str(my_node) + ".psms")
    current_working_psm_set.write_output(open(output_filename, "w"), True)






if __name__ == "__main__":
    main()
