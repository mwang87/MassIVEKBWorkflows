#!/usr/bin/python


import sys
import getopt
import os
import ming_protein_library
import ming_psm_library
import ming_parallel_library
import library_creation
import json

# Given first pass search results and second pass search results
# Also the second pass proteins

def usage():
    usage_string = "<first pass search results> <second pass search results> <fasta database> <second pass proteins> <Output Peptides First Pass> <Output Peptides Second Pass> <output PSMs list with updated evalues>"
    usage_string += " <output 5% FDR original psms> <output 5% FDR updated eval psms>"
    print(usage_string)

def map_sequence_to_proteins(input_list):
    return input_list[0].get_proteins_with_sequence(input_list[1])


def get_first_pass_variant_set(first_pass_results_filename):
    psm_list_first_pass = ming_psm_library.PSMset(first_pass_results_filename)
    psm_list_first_pass.load_MSGF_Plus_tsvfile(first_pass_results_filename)
    psm_list_first_pass.filter_to_fdr_by_length(0.05)
    print "First Pass PSMs: " + str(len(psm_list_first_pass))

    full_peptides_list_first_pass = library_creation.create_library_unique_peptides_filtered([psm_list_first_pass], filter_by_length=True)
    print "First Pass Variants: " + str(len(full_peptides_list_first_pass))

    return full_peptides_list_first_pass

def get_second_pass_psms(second_pass_results_filename):
    psm_list_second_pass = ming_psm_library.PSMset(second_pass_results_filename)
    psm_list_second_pass.load_MSGF_Plus_tsvfile(second_pass_results_filename)

    return psm_list_second_pass

def output_updated_variant_evalues(full_peptides_list_first_pass, psm_list_second_pass, output_file_handle=None):
    second_pass_filename_scan_to_psm_map = {}
    for psm in psm_list_second_pass.psms:
        key = psm.filename + ":" + str(psm.scan)
        second_pass_filename_scan_to_psm_map[key] = psm

    for variant in full_peptides_list_first_pass.peptide_list:
        output_line = ""

        best_psm = variant.get_best_psm()
        key_lookup = best_psm.filename + ":" + best_psm.scan
        if key_lookup in second_pass_filename_scan_to_psm_map:
            second_pass_psm = second_pass_filename_scan_to_psm_map[key_lookup]
            if best_psm.annotation == second_pass_psm.annotation:
                #print "SAME"
                if second_pass_psm.sorting_value() > variant.sorting_value():
                    score_improvement = second_pass_psm.sorting_value() - variant.sorting_value()
                    #print "SCORE IMPROVEMENT: " + str(score_improvement)

                    output_line = second_pass_psm.annotation + "\t" + str(second_pass_psm.sorting_value()) + "\t"
                    output_line += str(variant.is_decoy()) + "\t" + str(variant.fdr) + "\t" + second_pass_psm.filename + "\t" + str(second_pass_psm.scan)
                    output_line += "\t" + str(score_improvement)
                else:
                    #print "SCORE SAME OR LESS"
                    output_line = str(variant) + "\t0"
            else:
                #print "NOT SAME"
                output_line = str(variant) + "\t0"
        else:
            output_line = str(variant) + "\t0"

        if output_file_handle == None:
            continue
        else:
            output_file_handle.write(output_line + "\n")

def update_evalues_first_second_pass(first_pass_psms, second_pass_psms):
    second_pass_psms_dict = {}

    #Creating the dicts
    for psm in second_pass_psms.psms:
        key = psm.filename + ":" + str(psm.scan) + ":" + psm.annotation
        second_pass_psms_dict[key] = psm

    print("Score improvement for PSMs")
    for psm in first_pass_psms.psms:
        psm.extra_metadata["score_improvement"] = str(0.0)

        key = psm.filename + ":" + str(psm.scan) + ":" + psm.annotation

        if key in second_pass_psms_dict:
            second_pass_score = second_pass_psms_dict[key].score
            second_pass_annotation = second_pass_psms_dict[key].annotation

            first_pass_annotation = psm.annotation
            first_pass_score = psm.score

            if first_pass_annotation == second_pass_annotation:
                if second_pass_score > first_pass_score:
                    psm.score = second_pass_score
                    psm.extra_metadata["score_improvement"] = str(second_pass_score - first_pass_score)



def update_psm_set_with_second_pass_psms(first_pass_psms, second_pass_psms, output_psms, FDR=0.05):
    #print(second_pass_psms)

    print("Loading second pass PSMs", second_pass_psms)
    psm_list_second_pass = ming_psm_library.PSMset(second_pass_psms)
    psm_list_second_pass.load_MSGF_Plus_tsvfile(second_pass_psms)
    psm_list_second_pass.remove_duplicated_rows()

    print("Loading first pass PSMs", first_pass_psms)
    psm_list_first_pass = ming_psm_library.PSMset(first_pass_psms)
    psm_list_first_pass.load_MSGF_Plus_tsvfile(first_pass_psms)
    psm_list_first_pass.remove_duplicated_rows()

    update_evalues_first_second_pass(psm_list_first_pass, psm_list_second_pass)

    psm_list_first_pass.filter_to_fdr_by_length(FDR)

    #Writing out the results
    psm_list_first_pass.write_output(open(output_psms, "w"), write_extra_metadata=True)


def main():
    first_pass_results_filename = sys.argv[1]
    second_pass_results_filename = sys.argv[2]
    fasta_db_filename = sys.argv[3]

    second_pass_proteins_filename = sys.argv[4]

    output_first_pass_peptides = sys.argv[5]
    output_second_pass_peptides = sys.argv[6]

    output_psms_first_pass = sys.argv[7]
    output_psms_updated_evalues = sys.argv[8]

    output_original_high_FDR_psms = sys.argv[9]
    output_updated_high_FDR_psms = sys.argv[10]

    #Low FDR Original and updated evals
    psm_list_first_pass = ming_psm_library.PSMset(first_pass_results_filename)
    psm_list_first_pass.load_MSGF_Plus_tsvfile(first_pass_results_filename)
    psm_list_first_pass.remove_duplicated_rows()
    psm_list_first_pass.filter_to_fdr_by_length(0.01)
    psm_list_first_pass.write_output(open(output_psms_first_pass, "w"))

    update_psm_set_with_second_pass_psms(first_pass_results_filename, second_pass_results_filename, output_psms_updated_evalues)


    #High FDR for other purposes to show things
    psm_list_first_pass = ming_psm_library.PSMset(first_pass_results_filename)
    psm_list_first_pass.load_MSGF_Plus_tsvfile(first_pass_results_filename)
    psm_list_first_pass.remove_duplicated_rows()
    psm_list_first_pass.filter_to_fdr_by_length(0.01)
    psm_list_first_pass.write_output(open(output_original_high_FDR_psms, "w"))

    update_psm_set_with_second_pass_psms(first_pass_results_filename, second_pass_results_filename, output_updated_high_FDR_psms, 0.05)


    #Precursor Level
    psm_list_first_pass = ming_psm_library.PSMset(first_pass_results_filename)
    psm_list_first_pass.load_MSGF_Plus_tsvfile(first_pass_results_filename)
    psm_list_first_pass.remove_duplicated_rows()
    full_peptides_list_first_pass = library_creation.create_library_unique_peptides_filtered([psm_list_first_pass], filter_by_length=True)
    full_peptides_list_first_pass.write_output(open(output_first_pass_peptides, "w"))

    psm_list_second_pass = ming_psm_library.PSMset(second_pass_results_filename)
    psm_list_second_pass.load_MSGF_Plus_tsvfile(second_pass_results_filename)
    psm_list_second_pass.remove_duplicated_rows()
    full_peptides_list_second_pass = library_creation.create_library_unique_peptides_filtered([psm_list_second_pass], filter_by_length=True)
    full_peptides_list_second_pass.write_output(open(output_second_pass_peptides, "w"))




if __name__ == "__main__":
    main()
