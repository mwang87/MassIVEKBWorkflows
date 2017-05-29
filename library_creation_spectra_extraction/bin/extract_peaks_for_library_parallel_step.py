#!/usr/bin/python


import sys
import getopt
import os
import json
import ming_psm_library
import library_creation
import ming_fileio_library
import ming_spectrum_library
import ming_proteosafe_library
import ming_numerical_utilities
import hashlib
import shutil
from collections import defaultdict

def usage():
    print("<paramxml> <input tsv file> <output folder>")

def group_psms_by_filename(psm_set):
    filename_to_psm_dict = defaultdict(list)
    for psm in psm_set.psms:
        filename = psm.filename
        filename_to_psm_dict[filename].append(psm)

    return filename_to_psm_dict

def extract_psms_from_filename(filename, psms_list, snr_threshold, minimum_explained_intensity, min_signal_peaks, min_number_of_peaks_within_1_percent_of_max, min_number_of_annotated_ions, max_ppm_error):
    full_path = os.path.join(ming_proteosafe_library.PROTEOSAFE_USER_UPLOADS_DIR, filename,)
    print("loading ", full_path)
    spectrum_collection = ming_spectrum_library.SpectrumCollection(full_path)
    spectrum_collection.load_from_file(drop_ms1=True)
    spectrum_list = []
    for psm in psms_list:
        scan = psm.scan
        protein = "PROTEIN"
        if psm.decoy == 1:
            protein = "CREATION_FALSE_PROTEIN"

        loaded_spectrum = spectrum_collection.scandict[scan]
        loaded_spectrum.filter_precursor_peaks()
        number_of_signal_peaks = loaded_spectrum.get_number_of_signal_peaks(SNR_Threshold=3)
        number_of_peaks_within_1_percent_of_max = loaded_spectrum.get_number_of_peaks_within_percent_of_max(percent=1.0)
        number_of_peaks_within_5_percent_of_max = loaded_spectrum.get_number_of_peaks_within_percent_of_max(percent=5.0)
        annotated_peak_count = ming_psm_library.calculated_number_annotated_peaks(loaded_spectrum.peaks, loaded_spectrum.charge, psm.get_annotation_without_charge(), 0.1)
        explained_intensity = ming_psm_library.calculated_explained_intensity(loaded_spectrum.peaks, loaded_spectrum.charge, psm.get_annotation_without_charge(), 0.1)
        number_of_ions_annotated_above_SNR = ming_spectrum_library.calculated_number_unique_ions_annotated_in_signal(loaded_spectrum.peaks, min(loaded_spectrum.charge, 3), psm.get_annotation_without_charge(), 0.1, SNR=3.0)

        theoretical_mz = ming_psm_library.calculate_theoretical_peptide_mass(psm.get_annotation_without_charge(), psm.charge)
        mass_difference = abs(theoretical_mz - loaded_spectrum.mz)
        ppm_error = (mass_difference / theoretical_mz) * 1000000
        parent_mass_error = mass_difference * psm.charge

        if snr_threshold > 0.9:
            loaded_spectrum.filter_noise_peaks(snr_threshold)

        output_spectrum_dict = {}
        output_spectrum_dict["filename"] = filename
        output_spectrum_dict["protein"] = protein
        output_spectrum_dict["scan"] = loaded_spectrum.scan
        output_spectrum_dict["peaks"] = json.dumps(loaded_spectrum.peaks)
        output_spectrum_dict["mz"] = loaded_spectrum.mz
        output_spectrum_dict["charge"] = psm.charge
        output_spectrum_dict["score"] = psm.score
        output_spectrum_dict["kl_score"] = float(psm.extra_metadata["kl_strict"])
        output_spectrum_dict["annotation"] = psm.get_annotation_without_charge()
        output_spectrum_dict["collision_energy"] = loaded_spectrum.collision_energy
        output_spectrum_dict["precursor_intensity"] = loaded_spectrum.precursor_intensity
        output_spectrum_dict["signal_peaks"] = number_of_signal_peaks
        output_spectrum_dict["number_of_peaks_within_1_percent_of_max"] = number_of_peaks_within_1_percent_of_max
        output_spectrum_dict["number_of_peaks_within_5_percent_of_max"] = number_of_peaks_within_5_percent_of_max
        output_spectrum_dict["annotated_peak_count"] = annotated_peak_count
        output_spectrum_dict["number_of_ions_annotated_above_SNR"] = number_of_ions_annotated_above_SNR
        output_spectrum_dict["explained_intensity"] = explained_intensity
        output_spectrum_dict["ppm_error"] = ppm_error
        output_spectrum_dict["parent_mass_error"] = parent_mass_error
        if "proteosafe_task" in psm.extra_metadata:
            output_spectrum_dict["proteosafe_task"] = psm.extra_metadata["proteosafe_task"]
        else:
            output_spectrum_dict["proteosafe_task"] = ""

        #TODO FILTER OUT SPECTRA HERE
        if output_spectrum_dict["signal_peaks"] < min_signal_peaks:
            continue
        if output_spectrum_dict["number_of_peaks_within_1_percent_of_max"] < min_number_of_peaks_within_1_percent_of_max:
            continue
        if output_spectrum_dict["explained_intensity"] < minimum_explained_intensity:
            continue
        if output_spectrum_dict["number_of_ions_annotated_above_SNR"] < min_number_of_annotated_ions:
            continue
        if ppm_error > max_ppm_error:
            continue

        spectrum_list.append(output_spectrum_dict)

    return spectrum_list

def get_snr_filter(param_obj):
    min_snr = 0.0
    try:
        min_snr = float(param_obj["min_snr_filter"][0])
    except KeyboardInterrupt:
        raise
    except:
        min_snr = 0.0

    return min_snr

def main():
    input_paramxml = sys.argv[1]
    input_tsv_filename = sys.argv[2]
    intermediate_output_folder = sys.argv[3]
    output_file_bins = int(sys.argv[4])

    params_obj = ming_proteosafe_library.parse_xml_file(open(input_paramxml))
    snr_threshold = get_snr_filter(params_obj)

    #Filtering Criteria
    minimum_explained_intensity = 0.0
    min_number_of_peaks_within_1_percent_of_max = 0.0
    min_signal_peaks = 0.0
    min_number_of_annotated_ions = 0.0
    max_kl_strict_score = 50
    max_ppm_error = 100000000

    try:
        minimum_explained_intensity = float(params_obj["min_explained_intensity"][0])
        min_number_of_peaks_within_1_percent_of_max = float(params_obj["min_number_of_peaks_within_1_percent_of_max"][0])
        min_signal_peaks = float(params_obj["min_signal_peaks"][0])
        min_number_of_annotated_ions = float(params_obj["min_number_of_annotated_ions"][0])
        max_kl_strict_score = float(params_obj["kl_strict_max"][0])
        if max_kl_strict_score == 0:
            max_kl_strict_score = 50
        max_ppm_error = float(params_obj["max_ppm_error"][0])
    except:
        print("exception")
        minimum_explained_intensity = 0.0
        min_number_of_peaks_within_1_percent_of_max = 0.0
        min_signal_peaks = 0.0
        max_kl_strict_score = 50


    #lets find the 1% variant point, and then the naive solution is to to take the top scoring one
    psm_set = ming_psm_library.PSMset("")
    psm_set.load_PSM_tsvfile(input_tsv_filename, load_extra_metadata=True)

    filename_to_psm_dict = group_psms_by_filename(psm_set)

    #All output files, we are going to bin them starting now
    output_filename_prefix = os.path.join(intermediate_output_folder, ming_fileio_library.get_filename_without_extension(os.path.basename(input_tsv_filename)) + "_partition_")
    output_files = {}
    output_files_number_spectra = {}
    for i in range(output_file_bins):
        output_filename = output_filename_prefix + str(i) + ".json"
        output_file = open(output_filename, "w")
        output_file.write("[")
        output_files[i] = output_file
        output_files_number_spectra[i] = 0


    for filename in filename_to_psm_dict:
        extracted_spectra = extract_psms_from_filename(filename, filename_to_psm_dict[filename], snr_threshold, minimum_explained_intensity, min_signal_peaks, min_number_of_peaks_within_1_percent_of_max, min_number_of_annotated_ions, max_ppm_error)
        for spectrum in extracted_spectra:
            hashed_index = int(hashlib.sha1(spectrum["annotation"].encode('utf-8')).hexdigest(), 16) % (output_file_bins)
            if output_files_number_spectra[hashed_index] == 0:
                output_files[hashed_index].write(json.dumps(spectrum) + "\n")
            else:
                output_files[hashed_index].write("," + json.dumps(spectrum) + "\n")
            output_files_number_spectra[hashed_index] += 1

    for i in range(output_file_bins):
        output_files[i].write("]")
        output_files[i].close()


if __name__ == "__main__":
    main()
