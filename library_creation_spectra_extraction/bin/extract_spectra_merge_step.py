#!/usr/bin/python


import sys
import getopt
import os
import ujson as json
#import json
import hashlib
import ming_fileio_library


def usage():
    print("<input folder> <output folder>")

def main():
    input_intermediate_folder = sys.argv[1]
    output_folder = sys.argv[2]
    number_of_bins = int(sys.argv[3])

    #Determining all output handles based on mass
    output_files = {}
    output_files_number_spectra = {}
    for i in range(number_of_bins):
        output_filename = os.path.join(output_folder, str(i) + ".json")
        output_file = open(output_filename, "w")
        output_file.write("[")
        output_files[i] = output_file
        output_files_number_spectra[i] = 0

    input_json_files = ming_fileio_library.list_files_in_dir(input_intermediate_folder)

    for json_filename in input_json_files:
        print("Loading", json_filename)
        spectrum_list = json.load(open(json_filename))
        for spectrum in spectrum_list:
            hashed_index = int(hashlib.sha1(spectrum["annotation"].encode('utf-8')).hexdigest(), 16) % (number_of_bins)
            if output_files_number_spectra[hashed_index] == 0:
                output_files[hashed_index].write(json.dumps(spectrum) + "\n")
            else:
                output_files[hashed_index].write("," + json.dumps(spectrum) + "\n")
            output_files_number_spectra[hashed_index] += 1

    for i in range(number_of_bins):
        output_files[i].write("]")
        output_files[i].close()


if __name__ == "__main__":
    main()
