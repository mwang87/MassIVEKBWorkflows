#!/usr/bin/python


import sys
import getopt
import os
import ming_fileio_library
from collections import defaultdict

def usage():
    print("<input file> <output file>")



def main():
    input_filename = sys.argv[1]
    output_filename = sys.argv[2]

    row_count, table_data = ming_fileio_library.parse_table_with_headers(input_filename)

    output_dict = defaultdict(list)

    max_fdr = 0.01


    for i in range(row_count):
        sequence = table_data["sequence"][i]
        modified_sequence = sequence[:-2]
        fdr = float(table_data["FDR"][i])
        if fdr > max_fdr:
            continue

        for key in table_data:
            output_dict[key].append(table_data[key][i])
        output_dict["modified_sequence"].append(modified_sequence)

    ming_fileio_library.write_dictionary_table_data(output_dict, output_filename)


if __name__ == "__main__":
    main()
