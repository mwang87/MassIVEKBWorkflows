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


def main():
    paramxml_filename = sys.argv[1]
    psms_input_file = sys.argv[2]
    kl_input_file = sys.argv[3]
    output_psms_file = sys.argv[4]

    parameters_obj = ming_proteosafe_library.parse_xml_file(open(paramxml_filename))


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


    #Since we don't support more fields in the psm object, we're going to read this file in again as a tsv file and add the columns as necessary
    psm_rows, psm_table_data = ming_fileio_library.parse_table_with_headers(psms_input_file)
    psm_table_data["kl_strict"] = []
    psm_table_data["kl_unstrict"] = []
    psm_table_data["kl_interpeak"] = []
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

    #Change C to C+57
    #if "cysteine_protease.cysteine" in parameters_obj:
    #    if parameters_obj["cysteine_protease.cysteine"][0] == "c57":
    #        #Lets replace all the cysteines
    #        for i in range(psm_rows):
    #            psm_table_data["sequence"][i] = psm_table_data["sequence"][i].replace("C", "C+57")


    ming_fileio_library.write_dictionary_table_data(psm_table_data, output_psms_file)


if __name__ == "__main__":
    main()
