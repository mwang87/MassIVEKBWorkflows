#!/usr/bin/python

import sys
import getopt
import os
import json
import ming_fileio_library
import ming_spectrum_library

def usage():
    print("<input file> <output merged tsv folder> <output merged mgf folder>")

def main():
    input_intermediate_file = sys.argv[1]
    output_tsv_folder = sys.argv[2]
    output_mgf_folder = sys.argv[3]
    output_sptxt_folder = sys.argv[4]

    all_input_files = ming_fileio_library.list_files_in_dir(input_intermediate_folder)

    library_spectrum_collection = ming_spectrum_library.SpectrumCollection("library spectra")

    all_json_spectra_list = json.load(open(input_intermediate_file))
    print("Loaded", input_intermediate_file, len(all_json_spectra_list))
    for library_spectrum in list_of_library_spectra:
        lib_spec = ming_spectrum_library.PeptideLibrarySpectrum("", 0, 0, library_spectrum["peaks"], library_spectrum["mz"], library_spectrum["charge"], library_spectrum["annotation"], library_spectrum["protein"])
        if "score" in library_spectrum:
            lib_spec.score = library_spectrum["score"]
        if "variant_score" in library_spectrum:
            lib_spec.variant_score = library_spectrum["variant_score"]
        if "spectra_to_consider" in library_spectrum:
            lib_spec.num_spectra = library_spectrum["spectra_to_consider"]
        if "ranking" in library_spectrum:
            lib_spec.spectrum_ranking = library_spectrum["ranking"]
        if "proteosafe_task" in library_spectrum:
            lib_spec.proteosafe_task = library_spectrum["proteosafe_task"]
        if "originalspectrum_filename" in library_spectrum:
            lib_spec.originalfile_filename = library_spectrum["originalspectrum_filename"]
        if "originalspectrum_scan" in library_spectrum:
            lib_spec.originalfile_scan = str(library_spectrum["originalspectrum_scan"])

        library_spectrum_collection.spectrum_list.append(lib_spec)

    output_mgf_filename = os.path.join(output_mgf_folder, os.path.splitext(os.path.basename(input_intermediate_file))[0] + ".mgf")
    output_tsv_filename = os.path.join(output_tsv_filename, os.path.splitext(os.path.basename(input_intermediate_file))[0] + ".tsv")
    output_sptxt_filename = os.path.join(output_tsv_filename, os.path.splitext(os.path.basename(input_intermediate_file))[0] + ".sptxt")

    library_spectrum_collection_split.save_to_mgf(open(output_mgf_filename, "w"))
    library_spectrum_collection_split.save_to_tsv(open(output_tsv_filename, "w"), output_mgf_filename)

    try:
        library_spectrum_collection.save_to_sptxt(open(output_sptxt_filename, "w"))
    except:
        traceback.print_exc(file=sys.stdout)
        print("MEH")
    





if __name__ == "__main__":
    main()
