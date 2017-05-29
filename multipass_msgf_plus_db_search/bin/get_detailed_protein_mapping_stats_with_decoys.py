#!/usr/bin/python


import sys
import getopt
import os
import ming_protein_library
import ming_psm_library
import library_creation
import json

def usage():
    print "<input fasta> <search results> <output proteins as list>"

def main():
    input_fasta_filename = sys.argv[1]
    input_searchresults_filename = sys.argv[2]
    output_proteins_as_list = sys.argv[3]

    proteome = ming_protein_library.parse_fasta_proteome_file(input_fasta_filename)

    #for protein in proteome.protein_list:
    #    print protein.protein

    psm_list = ming_psm_library.PSMset(input_searchresults_filename)
    psm_list.load_MSGF_Plus_tsvfile(input_searchresults_filename)

    full_peptides_list = library_creation.create_library_unique_peptides_filtered([psm_list], 0.01, filter_by_length=True)

    target_peptide_strings = []
    decoy_peptide_strings = []
    for peptide_obj in full_peptides_list.peptide_list:
        peptide_to_search = peptide_obj.get_stripped_sequence()
        if peptide_obj.is_decoy():
            decoy_peptide_strings.append(peptide_to_search[::-1])
        else:
            target_peptide_strings.append(peptide_to_search)

    protein_coverage_of_targets = proteome.get_proteins_with_number_of_peptides_covered_map(target_peptide_strings)
    protein_coverage_of_decoys = proteome.get_proteins_with_number_of_peptides_covered_map(decoy_peptide_strings)

    output_file = open(output_proteins_as_list, "w")
    output_file.write("protein\tdecoy_count\ttarget_count\ttotal_count\tlength\n")

    for protein in protein_coverage_of_targets:
        output_string = protein + "\t"
        output_string += str(protein_coverage_of_decoys[protein]) + "\t"
        output_string += str(protein_coverage_of_targets[protein]) + "\t"
        output_string += str(protein_coverage_of_targets[protein] + protein_coverage_of_decoys[protein]) + "\t"
        output_string += str(len(proteome.protein_map[protein].sequence)) + "\n"

        output_file.write(output_string)
    output_file.close()

if __name__ == "__main__":
    main()
