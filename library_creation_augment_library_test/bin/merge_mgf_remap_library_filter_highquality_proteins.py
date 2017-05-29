#!/usr/bin/python


import sys
import getopt
import os
import ming_fileio_library
import ming_spectrum_library

def usage():
    print("<input mgf folder> <input fdr and stats filename> <input peptide protein mapping file> <output mgf folder> <output tsv folder> <output filename prefix>")

def load_precursor_to_protein_mapping(input_filename):
    row_count, table_data = ming_fileio_library.parse_table_with_headers(input_filename)

    precursor_to_protein_map = {}
    for i in range(row_count):
        precursor_string = table_data["original_peptide"][i]
        protein_string = table_data["proteins_mapped"][i]
        precursor_to_protein_map[precursor_string] = protein_string

    return precursor_to_protein_map


def proteins_to_include(input_filename):
    row_count, table_data = ming_fileio_library.parse_table_with_headers(input_filename)

    protein_set = set()
    for i in range(row_count):
        fdr = float(table_data["fdr"][i])
        protein = table_data["protein"][i]
        number_of_non_overlapping_sequences = len(table_data["number_of_non_overlapping_sequences"][i])

        if fdr <= 0.01 and number_of_non_overlapping_sequences > 1:
            protein_set.add(protein)

    return protein_set

def main():
    input_folder = sys.argv[1]

    input_protein_fdr_filename = sys.argv[2]
    input_peptide_protein_mapping_filename = sys.argv[3]

    precursor_to_protein_map = load_precursor_to_protein_mapping(input_peptide_protein_mapping_filename)
    included_proteins = proteins_to_include(input_protein_fdr_filename)

    output_mgf_folder = sys.argv[4]
    output_tsv_folder = sys.argv[5]
    output_sptxt_folder = sys.argv[6]

    output_filename_prefix = sys.argv[7]

    input_files = ming_fileio_library.list_files_in_dir(input_folder)

    all_library_spectra = []
    for input_filename in input_files:
        temp_spectra = ming_spectrum_library.load_mgf_peptide_library(input_filename)
        print("loaded ", len(temp_spectra), "from", input_filename)
        for spectrum in temp_spectra:
            peptide = spectrum.peptide
            protein = spectrum.protein

            #All Proteins
            new_proteins_set = set(precursor_to_protein_map[peptide].split(";"))
            if len(new_proteins_set.intersection(included_proteins)) == 0:
                continue

            if protein != "CREATION_FALSE_PROTEIN":
                spectrum.protein = precursor_to_protein_map[peptide]

            all_library_spectra.append(spectrum)

    library_spectrum_collection_split = ming_spectrum_library.SpectrumCollection("library spectra")
    library_spectrum_collection_split.spectrum_list = all_library_spectra

    output_tsv_filename = os.path.join(output_tsv_folder, output_filename_prefix + ".tsv")
    output_mgf_filename = os.path.join(output_mgf_folder, output_filename_prefix + ".mgf")
    output_sptxt_filename = os.path.join(output_mgf_folder, output_filename_prefix + ".sptxt")

    library_spectrum_collection_split.save_to_mgf(open(output_mgf_filename, "w"))
    library_spectrum_collection_split.save_to_tsv(open(output_tsv_filename, "w"), output_mgf_filename)

    library_spectrum_collection_split.save_to_sptxt(open(output_sptxt_filename, "w"))





if __name__ == "__main__":
    main()
