#!/usr/bin/python


import sys
import getopt
import os
import ming_proteosafe_library
import ming_fileio_library
import json
from collections import defaultdict

def usage():
    print("<param xml> <output provenance json> <task summary> <augment summary> <searhc files summary>")



def main():
    params_obj = ming_proteosafe_library.parse_xml_file(open(sys.argv[1]))

    augment_task_id = params_obj["task"][0]

    all_tasks_output_dict = defaultdict(list)
    all_augments_output_dict = defaultdict(list)
    all_spectrum_files_output_dict = defaultdict(list)

    search_task_to_augment = {}
    search_task_to_extraction = {}

    all_search_tasks = set()

    process_tree = True
    while process_tree:
        print("AUGMENT", augment_task_id, len(augment_task_id))
        augment_task_information = ming_proteosafe_library.get_task_information("proteomics2.ucsd.edu", augment_task_id)

        extract_task_id = ""
        previous_augment_task_id = ""

        for filename in augment_task_information["files"]:
            if filename.find("unfiltered_peptide_list") != -1:
                previous_augment_task_id = ming_fileio_library.get_root_folder(filename.replace(ming_fileio_library.get_root_folder(filename) + "/", ""))
            if filename.find("extracted_spectra_peptides_merged") != -1:
                extract_task_id = ming_fileio_library.get_root_folder(filename.replace(ming_fileio_library.get_root_folder(filename) + "/", ""))

        previous_augment_task_id = previous_augment_task_id.strip()
        if len(previous_augment_task_id) < 10:
            process_tree = False

        print(previous_augment_task_id, extract_task_id)

        all_augments_output_dict["augment_task"].append(augment_task_id)
        all_augments_output_dict["extract_task"].append(extract_task_id)
        all_augments_output_dict["precursor_count"].append(0)
        all_augments_output_dict["timestamp"].append(augment_task_information["createtime"])


        #Processing extract task_id
        extract_task_info = ming_proteosafe_library.get_task_information("proteomics2.ucsd.edu", extract_task_id)
        extract_task_parameters = ming_proteosafe_library.get_task_parameters("proteomics2.ucsd.edu", extract_task_id)


        tasks_to_extract = json.loads(extract_task_parameters["tasks_to_consolidate"][0])

        for task in tasks_to_extract:
            search_task_to_augment[task] = augment_task_id
            search_task_to_extraction[task] = extract_task_id

            all_tasks_output_dict["search_task_id"].append(task)
            all_tasks_output_dict["extract_task_id"].append(extract_task_id)
            all_tasks_output_dict["augment_task_id"].append(augment_task_id)

            all_search_tasks.add(task)

        print(extract_task_parameters["task_file"][0])
        path_to_task_file = os.path.join("/data/ccms-data/uploads" , extract_task_parameters["task_file"][0][2:-1])
        if os.path.isfile(path_to_task_file):
            print("SEARCH FILE", path_to_task_file)
            try:
                row_count, table_data = ming_fileio_library.parse_table_with_headers(path_to_task_file)
                print("Rows", row_count)
                for i in range(row_count):
                    search_task_id = table_data["TASKID"][i]
                    print(i, search_task_id)

                    search_task_to_augment[search_task_id] = augment_task_id
                    search_task_to_extraction[search_task_id] = extract_task_id

                    all_tasks_output_dict["search_task_id"].append(search_task_id)
                    all_tasks_output_dict["extract_task_id"].append(extract_task_id)
                    all_tasks_output_dict["augment_task_id"].append(augment_task_id)

                    all_search_tasks.add(search_task_id)
            except:
                raise
                continue


        augment_task_id = previous_augment_task_id

    print(len(all_search_tasks))

    for i in range(len(all_tasks_output_dict["search_task_id"])):
        search_task = all_tasks_output_dict["search_task_id"][i]
        try:
            print(search_task)
            task_information = ming_proteosafe_library.get_task_information("proteomics2.ucsd.edu", search_task)
            all_tasks_output_dict["search_description"].append(task_information["description"])
            for filename in task_information["files"]:
                if filename.find(".mzXML") != -1 or filename.find(".mzML") != -1:
                    all_spectrum_files_output_dict["spectrum_filename"].append(filename)
                    all_spectrum_files_output_dict["search_task"].append(search_task)
                    all_spectrum_files_output_dict["search_description"].append(task_information["description"])
        except KeyboardInterrupt:
            raise
        except:
            all_tasks_output_dict["search_description"].append("")
            print("error", search_task)
            continue

    provenace_structure = {}
    provenace_structure["search_task_to_augment"] = search_task_to_augment
    provenace_structure["search_task_to_extraction"] = search_task_to_extraction

    open(sys.argv[2], "w").write(json.dumps(provenace_structure, indent=4))

    ming_fileio_library.write_dictionary_table_data(all_tasks_output_dict, sys.argv[3])
    ming_fileio_library.write_dictionary_table_data(all_augments_output_dict, sys.argv[4])
    ming_fileio_library.write_dictionary_table_data(all_spectrum_files_output_dict, sys.argv[5])




if __name__ == "__main__":
    main()
