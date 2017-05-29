#!/usr/bin/python

import ming_psm_library
import library_creation_utilities

#Taking all PSMs and then making them peptides and filtering 1%
#Then merging the peptides and checking the error rates
def create_library_unique_peptides_combined(all_psm_list):
    full_peptide_set = ming_psm_library.PeptideVariantSet("Combined")

    for psm_set in all_psm_list:
        peptide_set = ming_psm_library.PeptideVariantSet("Test")
        peptide_set.add_psms_set(psm_set)
        peptide_set.filter_to_fdr(0.01)

        full_peptide_set.add_variant_set(peptide_set)

    return full_peptide_set



#Taking all peptides and creating the union and then filtering
#at a fixed 1% FDR
def create_library_unique_peptides_filtered(all_psm_list, fdr=0.01, filter_by_length=False):
    full_peptide_set = ming_psm_library.PeptideVariantSet("Combined")

    for psm_set in all_psm_list:
        peptide_set = ming_psm_library.PeptideVariantSet("Test")
        peptide_set.add_psms_set(psm_set)

        full_peptide_set.add_variant_set(peptide_set)

    if filter_by_length:
        full_peptide_set.filter_to_fdr_by_length(fdr)
    else:
        full_peptide_set.filter_to_fdr(fdr)
    return full_peptide_set


#Taking all peptides and creating union and then filtering at a fixed 1% FDR at
#each peptide length
def create_library_merged_psm_list_separate_fdr_peptide_length(psm_set, fdr=0.01):
    full_peptide_set = ming_psm_library.PeptideVariantSet("Combined")


    peptide_length_map = {}
    for psm in psm_set.psms:
        peptide_length = len(psm.get_stripped_sequence())
        if not peptide_length in peptide_length_map:
            peptide_length_map[peptide_length] = ming_psm_library.PSMset("length" + str(peptide_length))
        peptide_length_map[peptide_length].psms.append(psm)

    #Lets do FDR on each length
    for peptide_length in peptide_length_map:
        #print peptide_length_map[peptide_length]
        peptide_set = ming_psm_library.PeptideVariantSet("Test")
        peptide_set.add_psms_set(peptide_length_map[peptide_length])
        peptide_set.filter_to_fdr(fdr)

        full_peptide_set.add_variant_set(peptide_set)

    #print full_peptide_set.peptide_list
    return full_peptide_set



#Taking all peptides and only keeping around peptides that have a
#replicate count that is greater than K
def create_library_minimum_replicate_count(all_psm_list):
    full_peptide_set = ming_psm_library.PeptideVariantSet("Combined")

    for psm_set in all_psm_list:
        peptide_set = ming_psm_library.PeptideVariantSet("Test")
        peptide_set.add_psms_set(psm_set)

        full_peptide_set.add_variant_set(peptide_set)

    full_peptide_set_new = ming_psm_library.PeptideVariantSet("Combined")

    for variant in full_peptide_set.peptide_list:
        if variant.get_spectrum_count() > 1:
            full_peptide_set_new.add_variant(variant)

    return full_peptide_set_new


#Taking each variant and checking variant consistency of psms and spectra
def create_library_variant_consistency(all_psm_list):
    full_peptide_set = ming_psm_library.PeptideVariantSet("Combined")

    for psm_set in all_psm_list:
        peptide_set = ming_psm_library.PeptideVariantSet("Test")
        peptide_set.add_psms_set(psm_set)

        full_peptide_set.add_variant_set(peptide_set)

    for variant in full_peptide_set.peptide_list:
        library_creation_utilities.check_variant_psm_consistency(variant)

    return full_peptide_set
