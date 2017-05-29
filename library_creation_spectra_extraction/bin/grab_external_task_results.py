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
from collections import defaultdict

def usage():
    print("<input params> <output peptide results folder> <output psms folder at 1% FDR>")


def grab_all_results(task_list, output_peptide_directory, output_psm_directory, output_summary_filename):
    results_list = []
    for task in task_list:
        ret_value = grab_single_result(task, output_peptide_directory, output_psm_directory)
        results_list.append(ret_value)

    summary_dictionary = defaultdict(list)
    for result in results_list:
        for key in result.keys():
            summary_dictionary[key].append(result[key])

    ming_fileio_library.write_dictionary_table_data(summary_dictionary, output_summary_filename)


def grab_single_result(task_id, output_peptide_directory, output_psm_directory):
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
        return_dict = grab_results_from_MSGFDB(task_id, user, output_peptide_directory, output_psm_directory)
        return return_dict

    if len(path_to_secondpass_peptides_files_list) == 1:
        return_dict = grab_results_from_multipass(task_id, user, output_peptide_directory, output_psm_directory)
        return return_dict

#Grabbing both the rescored peptide results from multi pass as well, returns number of PSMs and peptides
def grab_results_from_multipass(task_id, user, output_peptide_directory, output_psm_directory):
    return_dict = {}
    return_dict["number_psms"] = 0
    return_dict["number_peptides"] = 0
    return_dict["task_id"] = task_id

    #Copying the psm files
    path_to_psm_files_list = ming_proteosafe_library.get_proteosafe_result_file_path(task_id, user, "updated_eval_psms_with_kl_with_ambiguity")
    if len(path_to_psm_files_list) == 1:
        output_psm_path = os.path.join(output_psm_directory, task_id + ".psms")
        path_to_param_file = ming_proteosafe_library.get_proteosafe_result_file_path(task_id, user, "params")[0]

        #path_to_merged_results = ming_proteosafe_library.get_proteosafe_backend_result_file_path(task_id, "mergedResult", "proteomics2")[0]
        print(ming_proteosafe_library.get_proteosafe_result_file_path(task_id, user, "mergedResult"))
        path_to_merged_results = ming_proteosafe_library.get_proteosafe_result_file_path(task_id, user, "mergedResult")[0]

        print(path_to_psm_files_list[0] + " to " + output_psm_path)
        #name_demangle_filenames(path_to_psm_files_list[0], output_psm_path, path_to_param_file, "filename", "filename")
        name_demangle_filenames_and_instrument_collision(path_to_psm_files_list[0], output_psm_path, path_to_param_file, path_to_merged_results, "filename", "filename")

        #Now lets generate the peptide list from the psm list
        psm_set = ming_psm_library.PSMset("task results")
        psm_set.load_PSM_tsvfile(output_psm_path)
        output_peptide_path = output_psm_path = os.path.join(output_peptide_directory, task_id + ".peptides")

        peptide_variant_set = library_creation.create_library_unique_peptides_filtered([psm_set], 0.01)
        peptide_variant_set.write_output(open(output_peptide_path, "w"))

        return_dict["number_psms"] = len(psm_set.psms)
        return_dict["number_peptides"] = len(peptide_variant_set.peptide_list)

    return return_dict


def grab_results_from_MSGFDB(task_id, user, output_peptide_directory, output_psm_directory):
    return_dict = {}
    return_dict["number_psms"] = 0
    return_dict["number_peptides"] = 0
    return_dict["task_id"] = task_id

    path_to_psm_files_list = ming_proteosafe_library.get_proteosafe_result_file_path(task_id, user, "merged_result_with_kl_with_ambiguity")
    if len(path_to_psm_files_list) == 1:
        output_psm_path = os.path.join(output_psm_directory, task_id + ".psms")
        path_to_param_file = ming_proteosafe_library.get_proteosafe_result_file_path(task_id, user, "params")[0]

        #path_to_merged_results = ming_proteosafe_library.get_proteosafe_backend_result_file_path(task_id, "mergedResult", "proteomics2")[0]
        print(ming_proteosafe_library.get_proteosafe_result_file_path(task_id, user, "mergedResult"))
        path_to_merged_results = ming_proteosafe_library.get_proteosafe_result_file_path(task_id, user, "mergedResult")[0]

        print(path_to_psm_files_list[0] + " to " + output_psm_path)
        #name_demangle_filenames(path_to_psm_files_list[0], output_psm_path, path_to_param_file, "filename", "filename")
        name_demangle_filenames_and_instrument_collision(path_to_psm_files_list[0], output_psm_path, path_to_param_file, path_to_merged_results, "filename", "filename")

        #Now lets generate the peptide list from the psm list
        psm_set = ming_psm_library.PSMset("task results")
        psm_set.load_PSM_tsvfile(output_psm_path)
        output_peptide_path = os.path.join(output_peptide_directory, task_id + ".peptides")

        peptide_variant_set = library_creation.create_library_unique_peptides_filtered([psm_set], 0.01)
        peptide_variant_set.write_output(open(output_peptide_path, "w"))

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
    row_count, table_data = ming_fileio_library.parse_table_with_headers(input_file)
    mangled_mapping = ming_proteosafe_library.get_mangled_file_mapping(ming_proteosafe_library.parse_xml_file(open(path_to_param)))

    if not "FragMethod" in table_data:
        print("Demangling", path_to_original_results, input_file)
        collision_mapping = get_scan_mapping_for_collision_method(path_to_original_results)

        #Adding collision column
        table_data["FragMethod"] = []
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
    print(path_to_original_results)

    scan_header = "Scan#"
    if not scan_header in table_data:
        scan_header = "ScanNum"

    for i in range(row_count):
        key = table_data["#SpecFile"][i] + "_" + table_data[scan_header][i]
        mapping_dict[key] = table_data["FragMethod"][i]
    return mapping_dict


def main():
    params_filename = sys.argv[1]
    output_peptide_folder = sys.argv[2]
    output_psm_folder = sys.argv[3]
    output_summary = sys.argv[4]
    params_dict = ming_proteosafe_library.parse_xml_file(open(params_filename))

    source_tasks_text = params_dict["tasks_to_consolidate"][0]

    if len(source_tasks_text) > 0:
        source_tasks_list = json.loads(source_tasks_text)
        grab_all_results(source_tasks_list, output_peptide_folder, output_psm_folder, output_summary)
    else:
        open(output_summary, "w").write("None")



if __name__ == "__main__":
    main()
