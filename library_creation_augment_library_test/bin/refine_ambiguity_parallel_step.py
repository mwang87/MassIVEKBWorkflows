#!/usr/bin/python

import sys
import getopt
import os
import json
import ming_fileio_library
import ming_psm_library
import ming_proteosafe_library
from collections import defaultdict

def usage():
    print("<input file of result> <path to library identifications> <path to param.xml> <path to cutoff scores>")

def process_ambiguity(psm_list, mangled_mapping, library_scans_to_identification, cutoff_dict):
    output_results_dict = defaultdict(list)

    psm_by_filename_scan = defaultdict(list)

    #Grouping by scan and filename
    for psm in psm_list:
        spectrum_key = psm.filename + ":" + psm.scan
        psm_by_filename_scan[spectrum_key].append(psm)

    for key in psm_by_filename_scan:
        library_filename = psm_by_filename_scan[key][0].filename
        scan = psm_by_filename_scan[key][0].scan

        library_spectrum_key = library_filename + ":" + scan
        library_identification_object = library_scans_to_identification[library_spectrum_key]
        proteosafe_task = library_identification_object["proteosafe_task"]

        observed_annotations = set()
        observed_stripped_annotations = set()
        observed_stripped_annotations_scores = defaultdict(lambda: -1000)
        sequence_to_variant_map = defaultdict(list)


        for psm in psm_by_filename_scan[key]:
            annotation = psm.annotation
            stripped_annotation = ming_psm_library.strip_sequence(ming_psm_library.remove_charges_from_annotation(annotation))

            peptide_length = str(len(stripped_annotation))
            cutoff_score = 100000
            if proteosafe_task in cutoff_dict:
                task_cutoffs = cutoff_dict[proteosafe_task]
                if peptide_length in task_cutoffs:
                    cutoff_score = task_cutoffs[peptide_length]

            if psm.score >= cutoff_score:
                observed_annotations.add(annotation)
                observed_stripped_annotations.add(stripped_annotation)
                observed_stripped_annotations_scores[stripped_annotation] = max(psm.score, observed_stripped_annotations_scores[stripped_annotation])
                sequence_to_variant_map[stripped_annotation].append(annotation)

                #print(annotation, psm.score, cutoff_score)

        output_list = ["ALLSTATUS", str(key), str(len(observed_stripped_annotations)), library_identification_object["peptide"], library_identification_object["charge"]]
        #print("\t".join(output_list))

        ambiguity_category = "N/A"
        library_peptide = library_identification_object["peptide"]
        library_peptide_stripped = ming_psm_library.strip_sequence(ming_psm_library.remove_charges_from_annotation(library_identification_object["peptide"]))
        library_charge = library_identification_object["charge"]
        observed_annotations = len(observed_stripped_annotations)
        library_filename = library_identification_object["filename"]
        library_scan = library_identification_object["spectrumscan"]
        alternative_peptide = "N/A"

        if len(observed_stripped_annotations) == 2:
            ambiguous_stripped_sequence_1 = list(observed_stripped_annotations)[0]
            ambiguous_stripped_sequence_2 = list(observed_stripped_annotations)[1]

            ambiguous_mod_sequence_1 = sequence_to_variant_map[ambiguous_stripped_sequence_1][0]
            ambiguous_mod_sequence_2 = sequence_to_variant_map[ambiguous_stripped_sequence_2][0]

            sequence1_score = observed_stripped_annotations_scores[ambiguous_stripped_sequence_1]
            sequence2_score = observed_stripped_annotations_scores[ambiguous_stripped_sequence_2]

            score_delta = abs(sequence1_score - sequence2_score)

            #ambiguity_category = ming_ambiguity_library.categorize_peptide_distance(ambiguous_mod_sequence_1, ambiguous_mod_sequence_2)
            ambiguity_category = "N/A"

            if library_peptide_stripped == ambiguous_stripped_sequence_1:
                alternative_peptide = ambiguous_mod_sequence_2
            else:
                alternative_peptide = ambiguous_mod_sequence_1

            #print(ambiguous_mod_sequence_1, ambiguous_mod_sequence_2, ambiguity_category)

            #output_list = ["TWOAMBIGUOUS", str(key), str(len(observed_stripped_annotations)), library_identification_object["peptide"], library_identification_object["charge"], library_identification_object["score"], proteosafe_task, str(observed_stripped_annotations), str(cutoff_score), str(observed_annotations), str(ambiguity_category)]

            #print("\t".join(output_list))

        output_results_dict["ambiguity_category"].append(ambiguity_category)
        output_results_dict["library_peptide"].append(library_peptide)
        output_results_dict["library_peptide_stripped"].append(library_peptide_stripped)
        output_results_dict["library_charge"].append(library_charge)
        output_results_dict["observed_annotations"].append(observed_annotations)
        output_results_dict["library_filename"].append(library_filename)
        output_results_dict["library_scan"].append(library_scan)
        output_results_dict["alternative_peptide"].append(alternative_peptide)

    return output_results_dict



def library_scans_to_identification_info(library_identifications_filename):
    library_scans_to_identification = {}

    row_count, table_data = ming_fileio_library.parse_table_with_headers(library_identifications_filename)
    for i in range(row_count):
        filename = table_data["mgf_filename"][i]
        scan = table_data["spectrumscan"][i]
        peptide = table_data["peptide"][i]
        charge = table_data["charge"][i]
        score = table_data["score"][i]
        task_id = table_data["proteosafe_task"][i]

        identification_dict = {}
        identification_dict["peptide"] = peptide
        identification_dict["charge"] = charge
        identification_dict["score"] = score
        identification_dict["proteosafe_task"] = task_id
        identification_dict["filename"] = filename
        identification_dict["spectrumscan"] = scan

        key = os.path.basename(filename) + ":" + scan
        library_scans_to_identification[key] = identification_dict

    return library_scans_to_identification


def main():
    input_file_of_tsv_results = sys.argv[1]
    input_params_xml_filename = sys.argv[2]
    input_library_identifications_filename = sys.argv[3]
    input_cutoff_scores = sys.argv[4]
    output_folder = sys.argv[5]

    output_filename = os.path.join(output_folder, os.path.basename(input_file_of_tsv_results))

    params_object = ming_proteosafe_library.parse_xml_file(open(input_params_xml_filename))
    mangled_mapping = ming_proteosafe_library.get_mangled_file_mapping(params_object)

    library_scans_to_identification = library_scans_to_identification_info(input_library_identifications_filename)

    cutoff_dict = json.loads(open(input_cutoff_scores).read())

    psm_list = ming_psm_library.parse_MSGFPlus_tsvfile(input_file_of_tsv_results)
    output_results_dict = process_ambiguity(psm_list, mangled_mapping, library_scans_to_identification, cutoff_dict)

    ming_fileio_library.write_dictionary_table_data(output_results_dict, output_filename)


if __name__ == "__main__":
    main()
