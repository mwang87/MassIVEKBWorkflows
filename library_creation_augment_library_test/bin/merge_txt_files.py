#!/usr/bin/python


import sys
import getopt
import os
import ming_fileio_library

def usage():
    print("<input txt folder> <output txt file>")

def main():
    input_folder = sys.argv[1]
    output_filename_folder = sys.argv[2]

    input_files = ming_fileio_library.list_files_in_dir(input_folder)

    extension = os.path.split(input_files[0])[1]
    output_filename = os.path.join(output_filename_folder, "merged" + extension)
    output_file = open(output_filename, "w")

    for input_file in input_files:
        for line in open(input_file):
            output_file.write(line)

        output_file.write("\n")

    output_file.close()



if __name__ == "__main__":
    main()
