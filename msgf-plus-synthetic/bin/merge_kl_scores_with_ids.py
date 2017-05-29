#!/usr/bin/python


import sys
import getopt
import os
import ming_psm_library
import ming_fileio_library
import ming_proteosafe_library
import json


def usage():
    print "<param.xml> <input psms> <input kl file> <output psms>"

"""Returns mangled names and lists of both sets"""
def determine_set_of_target_and_decoy_spectrum_files(parameters_obj):
    reverse_mangling_mapping = ming_proteosafe_library.get_reverse_mangled_file_mapping(parameters_obj)

    target_filename = parameters_obj["target_run_filename"][0]

    target_filename_list = []
    decoy_filename_list = []

    for filename in reverse_mangling_mapping:
        if filename == target_filename:
            target_filename_list.append(reverse_mangling_mapping[filename])
        else:
            decoy_filename_list.append(reverse_mangling_mapping[filename])

    return target_filename_list, decoy_filename_list

def main():
    print(sys.argv)
    paramxml_filename = sys.argv[1]
    psms_input_file = sys.argv[2]
    kl_input_file = sys.argv[3]

    output_psms_file = sys.argv[4]
    output_decoy_psms_file = sys.argv[5]

    parameters_obj = ming_proteosafe_library.parse_xml_file(open(paramxml_filename))

    target_filename_list, decoy_filename_list = determine_set_of_target_and_decoy_spectrum_files(parameters_obj)

    input_psm_set = ming_psm_library.PSMset("input psms")
    input_psm_set.load_PSM_tsvfile(psms_input_file, load_extra_metadata=True)

    decoy_psm_set = ming_psm_library.PSMset("decoy psms")
    decoy_psm_set.psms = input_psm_set.synthetic_psms_by_length_decoy_set(target_filename_list, decoy_filename_list)

    print("GETTING ALL SYNETHTIC with 0% FDR")
    input_psm_set.filter_synthetic_psms_by_length(target_filename_list, decoy_filename_list)



    row_count, kl_data = ming_fileio_library.parse_table_with_headers(kl_input_file)
    kl_dict = {}
    for i in range(row_count):
        filename = os.path.basename(kl_data["Filename"][i])
        scan = kl_data["Scan"][i]
        kl_strict = (kl_data["KL Strict"][i])
        kl_unstrict = (kl_data["KL"][i])
        interpeak_intensity = (kl_data["Interpeak intensity"][i])
        key = filename + ":" + str(scan)
        kl_dict[key] = {"kl_strict" : kl_strict, "kl_unstrict" : kl_unstrict, "kl_interpeak" : interpeak_intensity}

    output_file = open(output_psms_file, "w")
    input_psm_set.write_output(output_file, write_extra_metadata=True)
    decoy_psm_set.write_output(open(output_decoy_psms_file, "w"), write_extra_metadata=True)
    output_file.close()


    #Since we don't support more fields in the psm object, we're going to read this file in again as a tsv file and add the columns as necessary
    psm_rows, psm_table_data = ming_fileio_library.parse_table_with_headers(output_psms_file)

    psm_table_data["kl_strict"] = []
    psm_table_data["kl_unstrict"] = []
    psm_table_data["kl_interpeak"] = []

    psm_table_data["ambiguity_total_score"] = []
    psm_table_data["first_second_unique_ratio"] = []
    psm_table_data["first_unique_count"] = []
    psm_table_data["first_unique_intensity"] = []
    psm_table_data["numberpsms"] = []
    psm_table_data["second_unique_count"] = []
    psm_table_data["second_unique_intensity"] = []
    psm_table_data["spectrum_unique_key"] = []
    psm_table_data["modified_sequence"] = []



    for i in range(psm_rows):
        key = psm_table_data["filename"][i] + ":" + psm_table_data["scan"][i]
        if key in kl_dict:
            psm_table_data["kl_strict"].append(kl_dict[key]["kl_strict"])
            psm_table_data["kl_unstrict"].append(kl_dict[key]["kl_unstrict"])
            psm_table_data["kl_interpeak"].append(kl_dict[key]["kl_interpeak"])
        else:
            psm_table_data["kl_strict"].append(-1)
            psm_table_data["kl_unstrict"].append(-1)
            psm_table_data["kl_interpeak"].append(-1)

        #writing the ambiguity stuff, but just assuming no ambiguity
        psm_table_data["ambiguity_total_score"].append("-1")
        psm_table_data["first_second_unique_ratio"].append("-1")
        psm_table_data["first_unique_count"].append("-1")
        psm_table_data["first_unique_intensity"].append("-1")
        psm_table_data["numberpsms"].append(1)
        psm_table_data["second_unique_count"].append("-1")
        psm_table_data["second_unique_intensity"].append("-1")
        psm_table_data["spectrum_unique_key"].append(key)
        psm_table_data["modified_sequence"].append(psm_table_data["sequence"][i][:-2])



    ming_fileio_library.write_dictionary_table_data(psm_table_data, output_psms_file)

if __name__ == "__main__":
    main()
