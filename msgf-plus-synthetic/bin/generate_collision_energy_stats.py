#!/usr/bin/python


import sys
import getopt
import os
import ming_psm_library
import ming_fileio_library
import ming_proteosafe_library
import ming_spectrum_library
import json
from collections import defaultdict


def get_file_stats(input_filename):
    output_list = []

    spectrum_collection = ming_spectrum_library.SpectrumCollection(input_filename)

    try:
        spectrum_collection.load_from_file(drop_ms1=True)
    except KeyboardInterrupt:
        raise
    except:
        print("Cannot load", input_filename)

    for spectrum in spectrum_collection.spectrum_list:
        output_dict = {}
        output_dict["fragmentation"] = spectrum.fragmenation_method
        output_dict["filename"] = os.path.basename(input_filename)
        output_dict["collision_energy"] = spectrum.collision_energy
        output_dict["scan"] = spectrum.scan
        output_list.append(output_dict)

    return output_list


def usage():
    print "<param.xml> <spectrum file> <output json file>"


def main():
    file_stats = get_file_stats(sys.argv[2])
    open(sys.argv[3], "w").write(json.dumps(file_stats, indent=4))


if __name__ == "__main__":
    main()
