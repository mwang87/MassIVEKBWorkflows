#!/usr/bin/python


import sys
import getopt
import os
import json
import ming_proteosafe_library
import ming_fileio_library
import library_creation
import ming_psm_library
import shutil
import pickle
from collections import defaultdict

def usage():
    print("<parallel json> <input params> <task_id file> <output peptide results folder> <output psms folder at 1% FDR>")

def filter_psms_with_params(params_obj, psm_list):
    print("Filtering")

    new_psm_list = []
    max_kl_strict_score = 50

    try:
        max_kl_strict_score = float(params_obj["kl_strict_max"][0])
        if max_kl_strict_score == 0:
            max_kl_strict_score = 50
    except KeyboardInterrupt:
        raise
    except:
        max_kl_strict_score = 50

    #Remove Ambiguous Spectra
    max_ambiguity_score = 50
    try:
        if params_obj["FilterSpectraWithAmbiguityScores"][0] == "Yes":
            max_ambiguity_score = 0
    except KeyboardInterrupt:
        raise
    except:
        max_ambiguity_score = 50

    #Filter out instrument type
    filter_instrument_type = False
    instrument_type = ""
    try:
        if params_obj["instrument_collision"][0] == "HCD":
            filter_instrument_type = True
            instrument_type = "HCD"
        if params_obj["instrument_collision"][0] == "CID":
            filter_instrument_type = True
            instrument_type = "CID"
    except KeyboardInterrupt:
        raise
    except:
        filter_instrument_type = False

    seen_psm_keys = set()
    #Also filtering out duplicate rows because of proteins
    for psm in psm_list:
        if float(psm.extra_metadata["kl_strict"]) > max_kl_strict_score:
            continue

        if float(psm.extra_metadata["ambiguity_total_score"]) > max_ambiguity_score:
            continue

        if float(psm.charge) > 8:
            continue

        if filter_instrument_type == True:
            if psm.frag_method != instrument_type:
                continue

        psm_key = psm.filename + "." + str(psm.scan) + "." + psm.annotation
        if psm_key in seen_psm_keys:
            continue
        seen_psm_keys.add(psm_key)

        new_psm_list.append(psm)

    new_psm_list_top_scoring_per_spectrum = []
    top_scoring_psm_per_spectrum = {}
    for psm in new_psm_list:
        spectrum_key = psm.filename + "." + str(psm.scan)

        if not spectrum_key in top_scoring_psm_per_spectrum:
            top_scoring_psm_per_spectrum[spectrum_key] = psm

        if psm.score > top_scoring_psm_per_spectrum[spectrum_key].score:
            top_scoring_psm_per_spectrum[spectrum_key] = psm

    for spectrum_key in top_scoring_psm_per_spectrum:
        new_psm_list_top_scoring_per_spectrum.append(top_scoring_psm_per_spectrum[spectrum_key])


    output_psm_list = []

    #Filtering out double N-term modifications
    for psm in new_psm_list_top_scoring_per_spectrum:
        annotation = psm.annotation
        if annotation.find("+42.011-17.027") != -1:
            continue
        if annotation.find("+42.011Q+0.984") != -1:
            continue
        if annotation.find("+42.011N+0.984") != -1:
            continue
        if annotation.find("+43.006-17.027") != -1:
            continue

        output_psm_list.append(psm)

    return output_psm_list

#Takes a psm set makes a peptide set, then saves out the PSMs as a pickle for efficient loading
def save_psms_as_peptides(psm_set, output_peptide_path, fdr):
    peptide_variant_set = library_creation.create_library_unique_peptides_filtered([psm_set], fdr, filter_by_length=True)

    psm_set = ming_psm_library.PSMset("task results")
    for peptide in peptide_variant_set.peptide_list:
        psm_set.psms.append(peptide.get_best_psm())

    output_pickle = open(output_peptide_path, 'wb')
    pickle.dump(psm_set, output_pickle, pickle.HIGHEST_PROTOCOL)
    output_pickle.close()

    return peptide_variant_set

def grab_all_results(task_list, output_peptide_directory, output_psm_directory, output_summary_filename, params_obj):
    results_list = []
    for task in task_list:
        ret_value = grab_single_result(task, output_peptide_directory, output_psm_directory, params_obj)
        results_list.append(ret_value)

    summary_dictionary = defaultdict(list)
    for result in results_list:
        for key in result.keys():
            summary_dictionary[key].append(result[key])

    ming_fileio_library.write_dictionary_table_data(summary_dictionary, output_summary_filename)


def grab_single_result(task_id, output_peptide_directory, output_psm_directory, params_obj):
    return_dict = {}
    return_dict["number_psms"] = 0
    return_dict["number_peptides"] = 0
    return_dict["task_id"] = task_id

    task_info = ming_proteosafe_library.get_task_information("proteomics2.ucsd.edu", task_id)
    user = task_info["user"]
    if task_info["status"] == "FAILED":
        return return_dict

    #lets check whether whether this has the peptide output, if not we can create it
    path_to_secondpass_peptides_files_list = ming_proteosafe_library.get_proteosafe_result_file_path(task_id, user, "updated_eval_psms_with_kl_with_ambiguity")

    if len(path_to_secondpass_peptides_files_list) == 0:
        return_dict = grab_results_from_task(task_id, user, output_peptide_directory, output_psm_directory, params_obj, "merged_result_with_kl_with_ambiguity")
        #return_dict = grab_results_from_MSGFDB(task_id, user, output_peptide_directory, output_psm_directory, params_obj)
        return return_dict

    if len(path_to_secondpass_peptides_files_list) == 1:
        return_dict = grab_results_from_task(task_id, user, output_peptide_directory, output_psm_directory, params_obj, "updated_eval_psms_with_kl_with_ambiguity")
        return return_dict

#Grabbing both the rescored peptide results from multi pass as well, returns number of PSMs and peptides
def grab_results_from_task(task_id, user, output_peptide_directory, output_psm_directory, params_obj, folder_for_results):
    return_dict = {}
    return_dict["number_psms"] = 0
    return_dict["number_peptides"] = 0
    return_dict["task_id"] = task_id

    #Copying the psm files
    path_to_psm_files_list = ming_proteosafe_library.get_proteosafe_result_file_path(task_id, user, folder_for_results)
    if len(path_to_psm_files_list) == 1:
        output_psm_path = os.path.join(output_psm_directory, task_id + ".psms")
        path_to_param_file = ming_proteosafe_library.get_proteosafe_result_file_path(task_id, user, "params")[0]

        #These are original results that are from MSGF+ that includes the fragmentation method
        print(task_id, user, ming_proteosafe_library.get_proteosafe_result_file_path(task_id, user, "mergedResult"))
        path_to_merged_results = ming_proteosafe_library.get_proteosafe_result_file_path(task_id, user, "mergedResult")[0]

        print(path_to_psm_files_list[0] + " to " + output_psm_path)
        #name_demangle_filenames(path_to_psm_files_list[0], output_psm_path, path_to_param_file, "filename", "filename")
        name_demangle_filenames_and_instrument_collision(path_to_psm_files_list[0], output_psm_path, path_to_param_file, path_to_merged_results, "filename", "filename")

        #Now lets generate the peptide list from the psm list
        psm_set = ming_psm_library.PSMset("task results")
        psm_set.load_PSM_tsvfile(output_psm_path, True)
        print("PSM Count", len(psm_set.psms))
        psm_set.psms = filter_psms_with_params(params_obj, psm_set.psms)
        #Setting the task of each psm
        for psm in psm_set.psms:
            psm.extra_metadata["proteosafe_task"] = task_id

        print("PSM Count Filtered", len(psm_set.psms))
        psm_set.filter_to_fdr_by_length(0.05)

        output_pickle = open(output_psm_path, 'wb')
        pickle.dump(psm_set, output_pickle, pickle.HIGHEST_PROTOCOL)
        output_pickle.close()

        output_peptide_path = output_psm_path = os.path.join(output_peptide_directory, task_id + ".peptides")

        peptide_variant_set = save_psms_as_peptides(psm_set, output_peptide_path, 0.05)

        return_dict["number_psms"] = len(psm_set.psms)
        return_dict["number_peptides"] = len(peptide_variant_set.peptide_list)

    return return_dict

def name_demangle_filenames(input_file, output_file, path_to_param, old_filename_header, new_filename_header):
    row_count, table_data = ming_fileio_library.parse_table_with_headers(input_file)
    mangled_mapping = ming_proteosafe_library.get_mangled_file_mapping(ming_proteosafe_library.parse_xml_file(open(path_to_param)))

    if old_filename_header == new_filename_header:
        for i in range(row_count):
            mangled_name = table_data[old_filename_header][i]
            unmangled_name = mangled_mapping[mangled_name]
            table_data[new_filename_header][i] = unmangled_name
    else:
        table_data[new_filename_header] = []
        for i in range(row_count):
            mangled_name = table_data[old_filename_header][i]
            unmangled_name = mangled_mapping[mangled_name]
            table_data[new_filename_header].append(unmangled_name)


    ming_fileio_library.write_dictionary_table_data(table_data, output_file)

def name_demangle_filenames_and_instrument_collision(input_file, output_file, path_to_param, path_to_original_results, old_filename_header, new_filename_header):
    row_count, table_data = ming_fileio_library.parse_table_with_headers(input_file, skip_incomplete_lines=True)
    mangled_mapping = ming_proteosafe_library.get_mangled_file_mapping(ming_proteosafe_library.parse_xml_file(open(path_to_param)))

    if not "FragMethod" in table_data:
        print("Demangling", path_to_original_results, input_file)
        collision_mapping = get_scan_mapping_for_collision_method(path_to_original_results)

        #Adding collision column
        table_data["FragMethod"] = []
        print(len(table_data["filename"]), len(table_data["scan"]))
        for i in range(row_count):
            key = table_data["filename"][i] + "_" + table_data["scan"][i]
            if key in collision_mapping:
                table_data["FragMethod"].append(collision_mapping[key])
            else:
                table_data["FragMethod"].append("NO_COLLISION")

    if old_filename_header == new_filename_header:
        for i in range(row_count):
            mangled_name = table_data[old_filename_header][i]
            unmangled_name = mangled_mapping[mangled_name]
            table_data[new_filename_header][i] = unmangled_name
    else:
        table_data[new_filename_header] = []
        for i in range(row_count):
            mangled_name = table_data[old_filename_header][i]
            unmangled_name = mangled_mapping[mangled_name]
            table_data[new_filename_header].append(unmangled_name)

    ming_fileio_library.write_dictionary_table_data(table_data, output_file)

def get_scan_mapping_for_collision_method(path_to_original_results):
    mapping_dict = {}
    row_count, table_data = ming_fileio_library.parse_table_with_headers(path_to_original_results)
    print("original results", path_to_original_results)

    scan_header = "Scan#"
    if not scan_header in table_data:
        scan_header = "ScanNum"

    print(len(table_data["#SpecFile"]), len(table_data[scan_header]))

    for i in range(row_count):
        key = table_data["#SpecFile"][i] + "_" + table_data[scan_header][i]
        mapping_dict[key] = table_data["FragMethod"][i]
    return mapping_dict


def main():
    parallel_json = json.loads(open(sys.argv[1]).read())
    params_filename = sys.argv[2]
    task_id_file = sys.argv[3]
    output_peptide_folder = sys.argv[4]
    output_psm_folder = sys.argv[5]
    #output_summary = sys.argv[5]
    params_dict = ming_proteosafe_library.parse_xml_file(open(params_filename))

    source_tasks_text = params_dict["tasks_to_consolidate"][0]

    row_count, task_file_table = ming_fileio_library.parse_table_with_headers(task_id_file)

    my_node = parallel_json["node_partition"]
    total_node = parallel_json["total_paritions"]

    output_summary = os.path.join(sys.argv[6], str(my_node))

    if len(source_tasks_text) > 0:
        source_tasks_list = json.loads(source_tasks_text)
        source_tasks_list += task_file_table["TASKID"]
        source_tasks_list.sort()
        source_tasks_list = source_tasks_list[my_node::total_node]
        grab_all_results(source_tasks_list, output_peptide_folder, output_psm_folder, output_summary, params_dict)
    else:
        open(output_summary, "w").write("None")



if __name__ == "__main__":
    main()
