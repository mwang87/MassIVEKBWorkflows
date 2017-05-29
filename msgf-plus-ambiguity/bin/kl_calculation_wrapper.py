#!/usr/bin/python


import sys
import getopt
import os
import json
import ming_fileio_library
import ming_proteosafe_library
from collections import defaultdict
import subprocess

def usage():
    print "<input param.xml> <input folder for mzXML> <output KL files folder> <scratch folder> <path to KL binary executable> <path to isotopes table>"

def main():
    input_param = ming_proteosafe_library.parse_xml_file(open(sys.argv[1]))
    input_folder = sys.argv[2]
    output_file = sys.argv[3]
    scratch_folder = sys.argv[4]
    path_to_executable = sys.argv[5]
    path_to_isotopes_table = sys.argv[6]

    #parent_mass_tolerance = input_param[]
    parent_mass_tolerance = 0.05

    all_input_file_paths = ming_fileio_library.list_files_in_dir(input_folder)

    output_kl_intermediates = []
    for input_file in all_input_file_paths:
        output_kl_file = os.path.join(scratch_folder, os.path.basename(input_file) + ".kl")
        cmd = path_to_executable + " --input " + input_file + " --output_summary " + output_kl_file + " " + "--peak_tolerance " + str(parent_mass_tolerance) + " --isotope_file " + path_to_isotopes_table + "  >/dev/null 2>&1 " 
        print(cmd)
        os.system(cmd)
        #subprocess.call([cmd])
        output_kl_intermediates.append(output_kl_file)

    combined_table = defaultdict(list)
    for output_kl_file in output_kl_intermediates:
        row_count, table_data = ming_fileio_library.parse_table_with_headers(output_kl_file)
        for key in table_data:
            combined_table[key] += table_data[key]

    ming_fileio_library.write_dictionary_table_data(combined_table, output_file)





if __name__ == "__main__":
    main()
