#!/usr/bin/python


import sys
import getopt
import os
import json
import hashlib
import ming_fileio_library
from collections import defaultdict

def usage():
    print("<input folder> <output file>")

def main():
    input_intermediate_folder = sys.argv[1]
    output_file = sys.argv[2]

    output_dict = defaultdict(list)

    total_rows = 0
    input_filenames = ming_fileio_library.list_files_in_dir(input_intermediate_folder)
    for input_filename in input_filenames:
        if total_rows > 10000000:
            continue

        row_count, table_data = ming_fileio_library.parse_table_with_headers(input_filename)
        total_rows += row_count
        for i in range(row_count):
            for key in table_data:
                output_dict[key].append(table_data[key][i])

    ming_fileio_library.write_dictionary_table_data(output_dict, output_file)

if __name__ == "__main__":
    main()
