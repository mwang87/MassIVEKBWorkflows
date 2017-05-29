#!/usr/bin/python

import sys
import getopt
import os
import json
import hashlib
import ming_fileio_library
import ming_proteosafe_library
import spectrum_alignment
import ming_numerical_utilities
import ming_psm_library
import ming_spectrum_library
from collections import defaultdict
import time
import copy

def usage():
    output_usage = "<input param xml file> <parameter json for parallelism> <filtered peptide list> <length score cutoff dict> "
    output_usage += "<existing_library_by_peptide_bins> <extracted_peaks_to_peptide_bins> "

    #output
    output_usage += "<output library candidates folder>"

    print(output_usage)



#returns path to files
def determine_filenames_to_load(my_node_number, params_obj, path_to_existing_library, path_to_new_library_spectra):
    existing_library_filename = ""
    new_library_filename = ""

    basic_filename = str(my_node_number) + ".json"
    mangled_file_mapping = ming_proteosafe_library.get_mangled_file_mapping(params_obj)

    existing_library_files = ming_fileio_library.list_files_in_dir(path_to_existing_library)
    new_library_spectra_files = ming_fileio_library.list_files_in_dir(path_to_new_library_spectra)

    for filename in existing_library_files:
        base_filename = os.path.basename(filename)
        unmangled_name = mangled_file_mapping[base_filename]
        if os.path.basename(unmangled_name) == basic_filename:
            existing_library_filename = os.path.join(path_to_existing_library, base_filename)

    for filename in new_library_spectra_files:
        base_filename = os.path.basename(filename)
        unmangled_name = mangled_file_mapping[base_filename]
        if os.path.basename(unmangled_name) == basic_filename:
            new_library_filename = os.path.join(path_to_new_library_spectra, base_filename)



    return existing_library_filename, new_library_filename


def filter_out_spectra_to_top(all_spectra, total_spectra_to_keep, top_per_dataset):
    #Look to reduce redundancy
    already_observed_spectra = set()

    spectra_per_dataset = defaultdict(list)
    for spectrum in all_spectra:
        spectrum_unique_key = spectrum["filename"] + ":" + str(spectrum["scan"])
        if spectrum_unique_key in already_observed_spectra:
            continue
        else:
            already_observed_spectra.add(spectrum_unique_key)

        filename = spectrum["filename"]
        root_folder = ming_fileio_library.get_root_folder(filename)
        spectra_per_dataset[root_folder].append(spectrum)

    all_spectra = []
    for dataset in spectra_per_dataset:
        list_of_spectra = spectra_per_dataset[dataset]
        list_of_spectra = sorted(list_of_spectra, key=lambda spectrum: spectrum["score"], reverse=True)
        list_of_spectra = list_of_spectra[:top_per_dataset]
        all_spectra += list_of_spectra

    all_spectra = sorted(all_spectra, key=lambda spectrum: spectrum["score"], reverse=True)
    all_spectra = all_spectra[:total_spectra_to_keep]

    return all_spectra

def main():
    param_filename = sys.argv[1]
    choose_consensus_params_filename = sys.argv[2]
    existing_library_spectra_folder = sys.argv[3]
    new_library_spectra_folder = sys.argv[4]

    output_library_all_spectra_json_folder = sys.argv[5]


    #Deciding on how to create consensus
    params_obj = ming_proteosafe_library.parse_xml_file(open(param_filename))

    parallel_params = json.loads(open(choose_consensus_params_filename).read())
    total_node_count = parallel_params["total_paritions"]
    my_node_number = parallel_params["node_partition"]

    consensus_selection_method = params_obj["ConsensusChoice"][0]


    #determine filenames
    existing_library_filename, new_library_filename = determine_filenames_to_load(my_node_number, params_obj, existing_library_spectra_folder, new_library_spectra_folder)

    output_library_all_spectra_json_filename = os.path.join(output_library_all_spectra_json_folder, str(my_node_number) + ".json")
    output_library_all_spectra_json_file_handle = open(output_library_all_spectra_json_filename, "w")

    print("Writing", output_library_all_spectra_json_filename)

    print(len(existing_library_filename), existing_library_filename)
    print(len(new_library_filename), new_library_filename)

    library_spectra = []
    top_scoring_to_keep = 100
    top_per_dataset = 20

    #If we are starting from scratch, so no existing library
    if len(existing_library_filename) == 0 and len(new_library_filename) == 0:
        print("no files to load")
        exit(0)

    if len(existing_library_filename) == 0 and len(new_library_filename) != 0:
        print("New Only")
        input_spectrum_file_handle = open(new_library_filename)
        for line in input_spectrum_file_handle:
            all_spectra = json.loads(line)

            if len(all_spectra) == 0:
                continue

            #Filter to only top K scoring psms
            #all_spectra = sorted(all_spectra, key=lambda spectrum: spectrum["score"], reverse=True)
            #all_spectra = all_spectra[:top_scoring_to_keep]
            """Filtering intelligently"""
            all_spectra = filter_out_spectra_to_top(all_spectra, top_scoring_to_keep, top_per_dataset)

            annotation = all_spectra[0]["annotation"] + "." + str(all_spectra[0]["charge"])

            output_library_all_spectra_json_file_handle.write(json.dumps(all_spectra))
            output_library_all_spectra_json_file_handle.write("\n")

    if len(existing_library_filename) != 0 and len(new_library_filename) != 0:
        print("New and Old")

        #load both files and iterate through them
        new_library_file_handle = open(new_library_filename)
        existing_library_file_handle = open(existing_library_filename)

        new_library_current_spectra_string = new_library_file_handle.readline()
        existing_library_current_spectra_string = existing_library_file_handle.readline()

        new_library_current_spectra = []
        existing_library_current_spectra = []

        new_library_precursor = ""
        existing_library_precursor = ""

        parse_new_spectra = True
        parse_existing_spectra = True

        new_spectra_ended = False
        existing_spectra_ended = False

        #new_library_current_spectra = json.loads(new_library_file_handle.readline())
        #existing_library_current_spectra = json.loads(existing_library_file_handle.readline())

        #new_library_precursor = new_library_current_spectra[0]["annotation"] + "." + str(new_library_current_spectra[0]["charge"])
        #existing_library_precursor = existing_library_current_spectra[0]["annotation"] + "." + str(existing_library_current_spectra[0]["charge"])

        #print(new_library_precursor, existing_library_precursor)

        while True:
            print(len(existing_library_current_spectra_string), len(new_library_current_spectra_string))

            if len(new_library_current_spectra_string) == 0:
                new_spectra_ended = True
                parse_new_spectra = False
            if len(existing_library_current_spectra_string) == 0:
                existing_spectra_ended = True
                parse_existing_spectra = False

            if existing_spectra_ended == True and new_spectra_ended == True:
                break

            if parse_new_spectra == True:
                new_library_current_spectra = json.loads(new_library_current_spectra_string)
                new_library_precursor = new_library_current_spectra[0]["annotation"] + "." + str(new_library_current_spectra[0]["charge"])
            if parse_existing_spectra == True:
                existing_library_current_spectra = json.loads(existing_library_current_spectra_string)
                existing_library_precursor = existing_library_current_spectra[0]["annotation"] + "." + str(existing_library_current_spectra[0]["charge"])


            if new_library_precursor == existing_library_precursor:
                print("FOUND BOTH")

                all_spectra = []
                all_spectra += new_library_current_spectra
                all_spectra += existing_library_current_spectra
                #Filter to only top K scoring psms
                #all_spectra = sorted(all_spectra, key=lambda spectrum: spectrum["score"], reverse=True)
                #all_spectra = all_spectra[:top_scoring_to_keep]
                """Filtering intelligently"""
                all_spectra = filter_out_spectra_to_top(all_spectra, top_scoring_to_keep, top_per_dataset)


                annotation = all_spectra[0]["annotation"] + "." + str(all_spectra[0]["charge"])
                print(annotation, len(all_spectra))

                #Get new spectra
                new_library_current_spectra_string = new_library_file_handle.readline()
                existing_library_current_spectra_string = existing_library_file_handle.readline()
                parse_new_spectra = True
                parse_existing_spectra = True

                output_library_all_spectra_json_file_handle.write(json.dumps(all_spectra))
                output_library_all_spectra_json_file_handle.write("\n")

            elif (new_library_precursor < existing_library_precursor and new_spectra_ended == False) or existing_spectra_ended == True:
                print("FOUND NEW", existing_spectra_ended, new_spectra_ended)

                all_spectra = new_library_current_spectra
                #Filter to only top K scoring psms
                #all_spectra = sorted(all_spectra, key=lambda spectrum: spectrum["score"], reverse=True)
                #all_spectra = all_spectra[:top_scoring_to_keep]
                """Filtering intelligently"""
                all_spectra = filter_out_spectra_to_top(all_spectra, top_scoring_to_keep, top_per_dataset)

                annotation = all_spectra[0]["annotation"] + "." + str(all_spectra[0]["charge"])
                print(annotation, len(all_spectra))

                #Get new spectra
                new_library_current_spectra_string = new_library_file_handle.readline()
                parse_new_spectra = True

                output_library_all_spectra_json_file_handle.write(json.dumps(all_spectra))
                output_library_all_spectra_json_file_handle.write("\n")


            elif (new_library_precursor > existing_library_precursor and existing_spectra_ended == False) or new_spectra_ended == True:
                print("FOUND EXISTING")

                all_spectra = existing_library_current_spectra
                #Filter to only top K scoring psms
                #all_spectra = sorted(all_spectra, key=lambda spectrum: spectrum["score"], reverse=True)
                #all_spectra = all_spectra[:top_scoring_to_keep]
                """Filtering intelligently"""
                all_spectra = filter_out_spectra_to_top(all_spectra, top_scoring_to_keep, top_per_dataset)

                annotation = all_spectra[0]["annotation"] + "." + str(all_spectra[0]["charge"])
                print(annotation, len(all_spectra))

                #Get new spectra
                existing_library_current_spectra_string = existing_library_file_handle.readline()
                parse_existing_spectra = True

                output_library_all_spectra_json_file_handle.write(json.dumps(all_spectra))
                output_library_all_spectra_json_file_handle.write("\n")

            else:
                print("Problem with Ordering")




if __name__ == "__main__":
    main()
