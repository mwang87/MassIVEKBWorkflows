#!/usr/bin/python


import sys
import getopt
import os
import ming_psm_library
import ming_fileio_library
import ming_proteosafe_library
import ming_spectrum_library
import spectrum_alignment
import json


def usage():
    print "<param.xml> <input spectrum file> <mergedResult> <ms2_metadata_folder> <filtered_psms_output_file>"


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

def load_collision_energy_mapping(input_folder):
    scan_maps = {}

    all_files = ming_fileio_library.list_files_in_dir(input_folder)
    for input_file in all_files:
        print(input_file)
        list_of_metadata = json.loads(open(input_file).read())
        for metadata in list_of_metadata:
            filename = metadata["filename"]
            scan = metadata["scan"]
            collision_energy = metadata["scan"]

            key = filename + ":" + str(scan)
            scan_maps[key] = metadata

    return scan_maps

def filter_psms_to_acceptable_metadata(input_psm_set, scan_maps, parameters_obj):
    acceptable_collision_energies = []
    for acceptable_energy_string in parameters_obj["acceptable_collision_energies"][0].split(";"):
        acceptable_collision_energies.append(float(acceptable_energy_string))

    new_psms_set = []
    for psm in input_psm_set.psms:
        filename = psm.filename
        scan = psm.scan
        key = filename + ":" + str(scan)
        scan_metadata = scan_maps[key]
        collision_energy = scan_metadata["collision_energy"]
        if collision_energy in acceptable_collision_energies:
            psm.collision_energy = collision_energy
            new_psms_set.append(psm)

    input_psm_set.psms = new_psms_set

def get_psms_to_current_file(input_psm_set, current_filename):
    new_psms_list = []
    for psm in input_psm_set.psms:
        if psm.filename == os.path.basename(current_filename):
            new_psms_list.append(psm)

    return new_psms_list

def get_psms_to_target_file(input_psm_set, target_filename_list):
    new_psms_list = []
    for psm in input_psm_set.psms:
        if psm.filename in target_filename_list:
            new_psms_list.append(psm)

    return new_psms_list

def filtering_out_blacklisted_decoys(input_decoy_psms, list_of_peptides):
    output_decoys_list = []

    for psm in input_decoy_psms:
        stripped_sequence = psm.get_stripped_sequence()
        if not stripped_sequence in list_of_peptides:
            output_decoys_list.append(psm)

    print("Filtered blacklist from ", len(input_decoy_psms), "to", len(output_decoys_list))

    return output_decoys_list

def filtering_out_high_scoring_decoys(input_decoy_psms, input_target_psms, target_filename, other_filename):
    input_decoy_psms = sorted(input_decoy_psms, key=lambda psm: psm.sorting_value(), reverse=True)

    print(target_filename, other_filename)

    target_collection = ming_spectrum_library.SpectrumCollection(target_filename)
    target_collection.load_from_mzXML(drop_ms1=True)
    decoy_collection = ming_spectrum_library.SpectrumCollection(other_filename)
    decoy_collection.load_from_mzXML(drop_ms1=True)

    top_scoring_precursor_target = {}
    for psm in input_target_psms:
        annotation = psm.annotation
        charge = psm.charge
        key = annotation + "." + str(charge)
        top_psm = psm
        if key in top_scoring_precursor_target:
            top_psm = top_scoring_precursor_target[key]
        else:
            top_scoring_precursor_target[key] = psm

        if psm.score > top_psm.score:
            top_scoring_precursor_target[key] = psm

    output_decoys_list = []
    for psm in input_decoy_psms[:200]:
        annotation = psm.annotation
        charge = psm.charge
        key = annotation + "." + str(charge)
        if key in top_scoring_precursor_target:
            print(key, psm.score, psm.scan)
            print(annotation, annotation[:-2])
            decoy_spectrum = decoy_collection.scandict[int(psm.scan)]
            target_spectrum = target_collection.scandict[int(top_scoring_precursor_target[key].scan)]
            cosine_score = spectrum_alignment.score_alignment_annotated_ion_peaks(decoy_spectrum.peaks, target_spectrum.peaks, 0, 0, 0.1, annotation, annotation, min(3, charge))
            print(cosine_score)
            if cosine_score < 0.7:
                output_decoys_list.append(psm)
        else:
            output_decoys_list.append(psm)

    for psm in input_decoy_psms[200:]:
        output_decoys_list.append(psm)

    return output_decoys_list

def filtering_redundant_identifications_per_scan(input_psm_list):
    output_psm_list = []
    set_already_included_spectrum_key = set()

    for psm in input_psm_list:
        spectrum_key = psm.filename + ":" + str(psm.scan)
        if spectrum_key in set_already_included_spectrum_key:
            continue
        else:
            set_already_included_spectrum_key.add(spectrum_key)
            output_psm_list.append(psm)

    return output_psm_list



def main():
    paramxml_filename = sys.argv[1]
    input_spectrum_filename = sys.argv[2]
    input_spectrum_all = sys.argv[3]
    psms_input_file = sys.argv[4]
    input_collision_energy_folder = sys.argv[5]
    output_psms_file = sys.argv[6]

    parameters_obj = ming_proteosafe_library.parse_xml_file(open(paramxml_filename))
    scan_metadata_maps = load_collision_energy_mapping(input_collision_energy_folder)

    target_filename_list, decoy_filename_list = determine_set_of_target_and_decoy_spectrum_files(parameters_obj)

    input_psm_set = ming_psm_library.PSMset("input psms")
    input_psm_set.load_MSGF_Plus_tsvfile(psms_input_file)

    """Filtering on Collision Energy"""
    print("Size Before Filtering", len(input_psm_set.psms))
    filter_psms_to_acceptable_metadata(input_psm_set, scan_metadata_maps, parameters_obj)
    print("Size After CE Filtering", len(input_psm_set.psms))

    """Filtering to current file"""
    current_file_psms = get_psms_to_current_file(input_psm_set, input_spectrum_filename)
    target_file_psms = get_psms_to_target_file(input_psm_set, target_filename_list)
    print(len(current_file_psms), len(target_file_psms))

    output_decoys_list = []
    if os.path.basename(input_spectrum_filename) in target_filename_list:
        #no filtering, just save
        print("Target")
        output_decoys_list = target_file_psms
    else:
        #Find top scoring hit for each precursor

        blacklisted_decoy_peptides = json.loads(parameters_obj["blacklisted_decoy_peptides_json"][0])
        current_file_psms = filtering_out_blacklisted_decoys(current_file_psms, blacklisted_decoy_peptides)
        output_decoys_list = filtering_out_high_scoring_decoys(current_file_psms, target_file_psms, os.path.join(input_spectrum_all, target_filename_list[0]), input_spectrum_filename)

    output_decoys_list = filtering_redundant_identifications_per_scan(output_decoys_list)
    input_psm_set.psms = output_decoys_list

    input_psm_set.write_output(open(output_psms_file, "w"))


if __name__ == "__main__":
    main()
