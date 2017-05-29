#!/usr/bin/python

import sys
import traceback
import getopt
import os
import json
import hashlib
import ming_fileio_library
import ming_proteosafe_library
import spectrum_alignment
import ming_numerical_utilities
import ming_psm_library
import ming_spectrum_library
from collections import defaultdict
import time
import copy

def usage():
    output_usage = "<input param xml file> <parameter json for parallelism> <filtered peptide list> <length score cutoff dict> "
    output_usage += "<merged peptide bins>"

    #output
    output_usage += "<all_library_spectrum_output_folder> <output library candidates folder>"

    print(output_usage)




#Takes a list of candidate spectra and summarizes it in output dict
def summarize_candidate_library_spectra(all_spectra, output_dict):
    for spectrum in all_spectra:
        most_similar_score = 0

        for other_spectrum in spectrum["score_dict"]:
            most_similar_score += spectrum["score_dict"][other_spectrum]
        if len(spectrum["score_dict"]) > 0:
            most_similar_score = most_similar_score/float(len(spectrum["score_dict"]))

        output_dict["filename"].append(spectrum["filename"])
        output_dict["scan"].append(spectrum["scan"])
        output_dict["mz"].append(spectrum["mz"])
        output_dict["score"].append(spectrum["score"])
        output_dict["charge"].append(spectrum["charge"])
        output_dict["annotation"].append(spectrum["annotation"])
        output_dict["number_of_ions_annotated_above_SNR"].append(spectrum["number_of_ions_annotated_above_SNR"])
        output_dict["annotated_peak_count"].append(spectrum["annotated_peak_count"])
        output_dict["ppm_error"].append(spectrum["ppm_error"])
        output_dict["explained_intensity"].append(spectrum["explained_intensity"])
        output_dict["most_similar_score"].append(most_similar_score)
        output_dict["kl_score"].append(spectrum["kl_score"])
        output_dict["number_of_peaks_within_1_percent_of_max"].append(spectrum["number_of_peaks_within_1_percent_of_max"])
        output_dict["precursor_intensity"].append(spectrum["precursor_intensity"])
        if "proteosafe_task" in spectrum:
            output_dict["proteosafe_task"].append(spectrum["proteosafe_task"])
        else:
            output_dict["proteosafe_task"].append("")



def choose_representative_spectrum_most_similary_combination_score(spectrum_list, minimum_spectra_for_combination=4):
    #Find the average spectrum, and then find the real spectrum that is closest to that average spectrum
    best_spectrum = None
    most_similar_score = -1000.0

    #for efficiency
    #spectrum_list = sorted(spectrum_list, key=lambda spectrum: spectrum["score"], reverse=True)
    total_scores_to_consider = len(spectrum_list)
    for spectrum in spectrum_list:
        spectrum_unique_key = spectrum["filename"] + ":" + str(spectrum["scan"])
        new_score_dict = {}
        existing_score_dict = {}

        if "score_dict" in spectrum:
            existing_score_dict = spectrum["score_dict"]

        all_scores = []
        for other_spectrum in spectrum_list:
            other_spectrum_unique_key = other_spectrum["filename"] + ":" + str(other_spectrum["scan"])

            if spectrum_unique_key == other_spectrum_unique_key:
                continue

            total_score = 0.0

            if other_spectrum_unique_key in existing_score_dict:
                total_score = existing_score_dict[other_spectrum_unique_key]
            else:
                total_score, reported_alignments = spectrum_alignment.score_alignment(spectrum["peaks"], other_spectrum["peaks"], spectrum["mz"], other_spectrum["mz"], 0.1)

            new_score_dict[other_spectrum_unique_key] = total_score
            all_scores.append(total_score)

        average_score = sum(all_scores) / float(total_scores_to_consider)
        if len(spectrum_list) < minimum_spectra_for_combination:
            explained_intensity = spectrum["explained_intensity"]
            annotated_ions = spectrum["number_of_ions_annotated_above_SNR"]
            average_score = average_score * explained_intensity * annotated_ions

        if average_score > most_similar_score:
            most_similar_score = average_score
            best_spectrum = spectrum

        #Saving the new score matrix
        spectrum["score_dict"] = new_score_dict

    #print(best_spectrum)
    return best_spectrum


#returns path to file to load as well as your parallel skipping function
def determine_filenames_to_load(my_node_number, total_parallel, path_to_merged_library_spectra):
    merged_library_filename = ""

    merged_library_files = ming_fileio_library.list_files_in_dir(path_to_merged_library_spectra)

    total_number_of_json_files = len(merged_library_files)

    json_file_number_to_load = my_node_number % total_number_of_json_files

    merged_library_filename = os.path.join(path_to_merged_library_spectra, str(json_file_number_to_load) + ".json")

    total_nodes_for_file = int(float(total_parallel)/float(total_number_of_json_files))
    if total_parallel % total_number_of_json_files > my_node_number % total_number_of_json_files:
        total_nodes_for_file += 1

    my_position_for_file = int(float(my_node_number)/float(total_number_of_json_files))

    return merged_library_filename, my_position_for_file, total_nodes_for_file

def create_library_spectrum(all_spectra, consensus_selection_method, score_cutoff_by_length, variant_to_score, library_candidates_output_dict, filter_peaks=False):
    representative_spectrum = None

    spectra_to_consider = []
    sequence = all_spectra[0]["annotation"]
    stripped_sequence = ming_psm_library.strip_sequence(sequence)
    length = len(stripped_sequence)
    score_cutoff = score_cutoff_by_length[length] - 0.01 #delta is for floating point errors
    for spectrum in all_spectra:
        if spectrum["score"] < score_cutoff:
            continue
        else:
            spectra_to_consider.append(spectrum)

    print("DEBUG", sequence, len(all_spectra), len(spectra_to_consider), score_cutoff)

    #Decode all the spectrum peaks
    for spectrum in spectra_to_consider:
        spectrum["peaks"] = json.loads(spectrum["peaks"])

    if consensus_selection_method == "MostSimilar_Combination_Score":
        representative_spectrum = choose_representative_spectrum_most_similary_combination_score(spectra_to_consider)

    #Summarizing
    summarize_candidate_library_spectra(spectra_to_consider, library_candidates_output_dict)

    representative_spectrum = copy.deepcopy(representative_spectrum)

    #Reencode peaks
    #for spectrum in all_spectra:
    #    spectrum["peaks"] = json.dumps(spectrum["peaks"])

    #Filtering out noise in library spectra
    if filter_peaks == True:
        representative_spectrum["peaks"] = ming_spectrum_library.filter_to_top_peaks(representative_spectrum["peaks"], 100)

    representative_ranking = 0
    representative_score = representative_spectrum["score"]
    for spectrum in spectra_to_consider:
        if spectrum["score"] >= representative_score:
            representative_ranking += 1

    #Creating library spectra
    library_spectrum = {}
    library_spectrum["peaks"] = representative_spectrum["peaks"]
    library_spectrum["charge"] = representative_spectrum["charge"]
    library_spectrum["annotation"] = representative_spectrum["annotation"]
    library_spectrum["mz"] = representative_spectrum["mz"]
    library_spectrum["protein"] = representative_spectrum["protein"]
    library_spectrum["score"] = representative_spectrum["score"]
    library_spectrum["spectra_to_consider"] = len(spectra_to_consider)
    library_spectrum["ranking"] = representative_ranking
    library_spectrum["originalspectrum_filename"] = representative_spectrum["filename"]
    library_spectrum["originalspectrum_scan"] = representative_spectrum["scan"]

    variant_key = representative_spectrum["annotation"] + "." + str(representative_spectrum["charge"])
    library_spectrum["variant_score"] = variant_to_score[variant_key]
    if "proteosafe_task" in representative_spectrum:
        library_spectrum["proteosafe_task"] = representative_spectrum["proteosafe_task"]
    else:
        library_spectrum["proteosafe_task"] = ""

    return library_spectrum

#Loads the annotation key, which is annotation dot charge as a set
def load_filtered_peptide_set(filtered_filename):
    row_count, table_data = ming_fileio_library.parse_table_with_headers(filtered_filename)
    print("Number of Variant Sequences", row_count)
    return set(table_data["variant_sequence"])

def load_score_cutoff_by_length(filtered_filename):
    row_count, table_data = ming_fileio_library.parse_table_with_headers(filtered_filename)

    score_cutoff_by_length = defaultdict(lambda: 10000)
    for i in range(row_count):
        length = int(table_data["length"][i])
        score = float(table_data["score"][i])
        score_cutoff_by_length[length] = min(score, score_cutoff_by_length[length])

    return score_cutoff_by_length

def load_variant_to_score(filtered_filename):
    row_count, table_data = ming_fileio_library.parse_table_with_headers(filtered_filename)

    variant_to_score = {}
    for i in range(row_count):
        variant = table_data["variant_sequence"][i]
        score = float(table_data["score"][i])
        variant_to_score[variant] = score

    return variant_to_score

def filter_out_spectra_to_top(all_spectra, total_spectra_to_keep, top_per_dataset):
    #Look to reduce redundancy
    already_observed_spectra = set()

    spectra_per_dataset = defaultdict(list)
    for spectrum in all_spectra:
        spectrum_unique_key = spectrum["filename"] + ":" + str(spectrum["scan"])
        if spectrum_unique_key in already_observed_spectra:
            continue
        else:
            already_observed_spectra.add(spectrum_unique_key)

        filename = spectrum["filename"]
        root_folder = ming_fileio_library.get_root_folder(filename)
        spectra_per_dataset[root_folder].append(spectrum)

    all_spectra = []
    for dataset in spectra_per_dataset:
        list_of_spectra = spectra_per_dataset[dataset]
        list_of_spectra = sorted(list_of_spectra, key=lambda spectrum: spectrum["score"], reverse=True)
        list_of_spectra = list_of_spectra[:top_per_dataset]
        all_spectra += list_of_spectra

    all_spectra = sorted(all_spectra, key=lambda spectrum: spectrum["score"], reverse=True)
    all_spectra = all_spectra[:total_spectra_to_keep]

    return all_spectra

def main():
    param_filename = sys.argv[1]
    choose_consensus_params_filename = sys.argv[2]
    filtered_peptide_list_filename = sys.argv[3]
    length_score_cutoff_filename = sys.argv[4]
    provenance_json_filename = sys.argv[5]
    merged_library_spectra_folder = sys.argv[6]

    output_library_json_folder = sys.argv[7]
    output_candidate_spectra_tsv_folder = sys.argv[8]

    filtered_peptide_set = load_filtered_peptide_set(filtered_peptide_list_filename)
    score_cutoff_by_length = load_score_cutoff_by_length(filtered_peptide_list_filename)
    variant_to_score = load_variant_to_score(filtered_peptide_list_filename)

    #Deciding on how to create consensus
    params_obj = ming_proteosafe_library.parse_xml_file(open(param_filename))

    parallel_params = json.loads(open(choose_consensus_params_filename).read())
    total_node_count = parallel_params["total_paritions"]
    my_node_number = parallel_params["node_partition"]

    consensus_selection_method = params_obj["ConsensusChoice"][0]

    #output dict for listing all candidates
    library_candidates_output_dict = defaultdict(list)

    #determine filenames
    merged_library_filename, my_position_for_file, total_nodes_for_file = determine_filenames_to_load(my_node_number, total_node_count, merged_library_spectra_folder)


    print(merged_library_filename, my_position_for_file, total_nodes_for_file)

    library_spectra = []

    input_spectrum_file_handle = open(merged_library_filename)
    line_count = 0
    for line in input_spectrum_file_handle:
        line_count += 1
        if line_count % total_nodes_for_file != my_position_for_file:
            #print("Should Skip")
            continue
        else:
            print("NOT SKIP")

        all_spectra = json.loads(line)

        if len(all_spectra) == 0:
            continue

        annotation = all_spectra[0]["annotation"] + "." + str(all_spectra[0]["charge"])
        print(annotation, len(all_spectra))

        if not annotation in filtered_peptide_set:
            continue

        library_spectrum = create_library_spectrum(all_spectra, consensus_selection_method, score_cutoff_by_length, variant_to_score, library_candidates_output_dict)
        library_spectra.append(library_spectrum)


    json.dump(library_spectra, open(os.path.join(output_library_json_folder, str(my_node_number) + ".json"), "w"))

    #Provenance Records
    provenance_records = json.loads(open(provenance_json_filename).read())

    #Modifying the output candidate file
    for i in range(len(library_candidates_output_dict["filename"])):
        proteosafe_task = library_candidates_output_dict["proteosafe_task"][i]
        if proteosafe_task in provenance_records["search_task_to_augment"]:
            library_candidates_output_dict["augment_task"].append(provenance_records["search_task_to_augment"][proteosafe_task])
        else:
            library_candidates_output_dict["augment_task"].append("")

        if proteosafe_task in provenance_records["search_task_to_extraction"]:
            library_candidates_output_dict["extract_task"].append(provenance_records["search_task_to_extraction"][proteosafe_task])
        else:
            library_candidates_output_dict["extract_task"].append("")

    #Outputting
    output_candidate_spectra_tsv_filename = os.path.join(output_candidate_spectra_tsv_folder, str(my_node_number) + ".tsv")
    ming_fileio_library.write_dictionary_table_data(library_candidates_output_dict, output_candidate_spectra_tsv_filename)


    """Converted Output"""
    output_tsv_folder = sys.argv[9]
    output_mgf_folder = sys.argv[10]
    output_sptxt_folder = sys.argv[11]

    library_spectrum_collection = ming_spectrum_library.SpectrumCollection("library spectra")

    for library_spectrum in library_spectra:
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

    output_mgf_filename = os.path.join(output_mgf_folder, str(my_node_number) + ".mgf")
    output_tsv_filename = os.path.join(output_tsv_folder, str(my_node_number) + ".tsv")
    output_sptxt_filename = os.path.join(output_sptxt_folder, str(my_node_number) + ".sptxt")

    library_spectrum_collection.save_to_mgf(open(output_mgf_filename, "w"))
    library_spectrum_collection.save_to_tsv(open(output_tsv_filename, "w"), output_mgf_filename)

    try:
        library_spectrum_collection.save_to_sptxt(open(output_sptxt_filename, "w"))
    except:
        traceback.print_exc(file=sys.stdout)
        print("MEH")





if __name__ == "__main__":
    main()
