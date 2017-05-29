#!/usr/bin/python


import sys
import getopt
import os
import ujson as json
#import json
import hashlib
import ming_fileio_library
from collections import defaultdict


def usage():
    print("<input json> <input folder> <output folder> <output peptide list folder>")

def main():
    input_json = json.loads(open(sys.argv[1]).read())
    input_intermediate_folder = sys.argv[2]
    output_folder = sys.argv[3]
    output_peptide_list_folder = sys.argv[4]

    my_node = input_json["node_partition"]

    output_filename = os.path.join(output_folder, str(my_node) + ".json")
    output_file = open(output_filename, "w")
    number_of_spectra = 0

    input_json_files = ming_fileio_library.list_files_in_dir(input_intermediate_folder)
    input_json_files.sort()

    all_spectra = []

    for json_filename in input_json_files:
        #Skip files
        json_basename = os.path.basename(json_filename).split(".")[0]
        bin_peptide = int(json_basename.split("_")[2])
        if bin_peptide != my_node:
            continue

        print("Loading", json_filename)
        spectrum_list = json.load(open(json_filename))
        all_spectra += spectrum_list
        print("Total Spectra", len(spectrum_list), len(all_spectra))

    peptide_dict = defaultdict(list)
    print("Creating hash")
    for spectrum in all_spectra:
        annotation = spectrum["annotation"] + "." + str(spectrum["charge"])
        peptide_dict[annotation].append(spectrum)

    print("writing out strings")
    all_annotation = list(peptide_dict.keys())
    all_annotation.sort()
    for annotation in all_annotation:
        output_file.write(json.dumps(peptide_dict[annotation]))
        output_file.write("\n")

    output_file.close()

    #Write out all the peptides into a file
    output_peptide_dict = defaultdict(list)
    for annotation_key in peptide_dict:
        max_score = -10
        if len(peptide_dict[annotation_key]) > 0:
            for spectrum in peptide_dict[annotation_key]:
                max_score = max(spectrum["score"], max_score)
            #max score per peptide
            output_peptide_dict["score"].append(max_score)
            output_peptide_dict["annotation_key"].append(annotation_key)
            output_peptide_dict["annotation"].append(peptide_dict[annotation_key][0]["annotation"])
            output_peptide_dict["charge"].append(peptide_dict[annotation_key][0]["charge"])
            output_peptide_dict["protein"].append(peptide_dict[annotation_key][0]["protein"])

    #writing out file
    output_peptide_filename = os.path.join(output_peptide_list_folder, str(my_node) + ".tsv")
    ming_fileio_library.write_dictionary_table_data(output_peptide_dict, output_peptide_filename)

if __name__ == "__main__":
    main()
