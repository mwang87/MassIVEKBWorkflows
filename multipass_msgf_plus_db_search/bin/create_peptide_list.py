#!/usr/bin/python


import sys
import getopt
import os
import ming_protein_library
import ming_psm_library
import library_creation
import json

def usage():
    print "<search results> <output peptide list> <output peptide list with decoy>"


def main():
    input_searchresults_filename = sys.argv[1]
    output_peptide_list = sys.argv[2]
    output_peptide_list_with_decoy_filename = sys.argv[3]

    psm_list = ming_psm_library.PSMset(input_searchresults_filename)
    psm_list.load_MSGF_Plus_tsvfile(input_searchresults_filename)
    psm_list.filter_to_fdr_by_length(0.01)
    print len(psm_list)

    full_peptides_list = library_creation.create_library_unique_peptides_filtered([psm_list], fdr=0.01, filter_by_length=True)

    output_file = open(output_peptide_list, "w")

    all_peptides = [peptide.get_stripped_sequence() for peptide in full_peptides_list.peptide_list]
    all_peptides = list(set(all_peptides))

    for peptide in all_peptides:
        output_file.write(peptide + "\n")


    #Now lets load the PSMs and keep all variants, and then output them with the decoys present
    print "GIVING US FULL RESULT SET"
    psm_list = ming_psm_library.PSMset(input_searchresults_filename)
    psm_list.load_MSGF_Plus_tsvfile(input_searchresults_filename)
    full_peptides_list = library_creation.create_library_unique_peptides_filtered([psm_list], 1.0)
    output_peptide_list_with_decoy_file = open(output_peptide_list_with_decoy_filename, "w")

    output_peptide_list_with_decoy_file.write(ming_psm_library.PeptideVariant.output_header() + "\n")
    for peptide in full_peptides_list.peptide_list:
        output_peptide_list_with_decoy_file.write(str(peptide) + "\n")








if __name__ == "__main__":
    main()
