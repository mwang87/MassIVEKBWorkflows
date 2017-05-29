#!/usr/bin/python


import sys
import getopt
import os
import ming_protein_library
import ming_psm_library
import library_creation
import json

def usage():
    print "<input fasta> <search results> <output proteins_files as json> <output proteins as list>"


def map_sequence_to_proteins(input_list):
    return input_list[0].get_proteins_with_sequence(input_list[1])

def main():
    input_fasta_filename = sys.argv[1]
    input_searchresults_filename = sys.argv[2]
    output_fasta_filename = sys.argv[3]
    output_proteins_as_list = sys.argv[4]

    proteome = ming_protein_library.parse_fasta_proteome_file(input_fasta_filename)

    #for protein in proteome.protein_list:
    #    print protein.protein

    psm_list = ming_psm_library.PSMset(input_searchresults_filename)
    psm_list.load_MSGF_Plus_tsvfile(input_searchresults_filename)
    psm_list.filter_to_fdr_by_length(0.01)
    print len(psm_list)

    full_peptides_list = library_creation.create_library_unique_peptides_filtered([psm_list], fdr=0.01, filter_by_length=True)

    #Testing efficient version fo this
    all_peptide_strings = []
    for peptide_obj in full_peptides_list.peptide_list:
        peptide_to_search = peptide_obj.get_stripped_sequence()
        all_peptide_strings.append(peptide_obj.get_stripped_sequence())

    all_proteins = proteome.get_proteins_covered_by_k_peptides(all_peptide_strings, 2, True)

    all_protein_names = []
    for protein in all_proteins:
        all_protein_names.append(protein.protein)

    output_protein_filename = output_fasta_filename
    open(output_protein_filename, "w").write(json.dumps(all_protein_names))

    #Outputting the list of proteins
    output_protein_list_file = open(output_proteins_as_list, "w")

    output_protein_list_file.write("Protein\n")
    for protein in all_protein_names:
        output_protein_list_file.write(protein + "\n")

    exit(0)




if __name__ == "__main__":
    main()
