#!/usr/bin/python


import sys
import getopt
import os
import ming_fileio_library
import json
from collections import defaultdict
import ming_spectrum_library
import ming_ambiguity_library
import ming_psm_library
import spectrum_alignment
import re

def usage():
    print("<input psms with kl> <input spectrum folder> <output file>")



def extract_annotated_peaks(ion_peak_mapping, peak_list, tolerance):
    extracted_peaks = []
    for peak in peak_list:
        mass = peak[0]
        for ion_peak in ion_peak_mapping:
            if abs(mass - ion_peak_mapping[ion_peak]) < tolerance:
                extracted_peaks.append(peak)
                break
    return extracted_peaks

def calculated_ambiguity(parameter_map, peak_tolerance):
    filename = parameter_map["filename"]
    scan_mapping = parameter_map["scan_mapping"]

    return_ambiguity_mapping = defaultdict(lambda: {})


    #Bypassing
    for scan in scan_mapping:
        score_summary = {}
        score_summary["ambiguity_total_score"] = -1
        score_summary["first_unique_count"] = -1
        score_summary["second_unique_count"] = -1
        score_summary["first_unique_intensity"] = -1
        score_summary["second_unique_intensity"] = -1
        score_summary["first_second_unique_ratio"] = -1

        return_ambiguity_mapping[scan] = score_summary

    return return_ambiguity_mapping

    spectrum_collection = ming_spectrum_library.SpectrumCollection(filename)
    spectrum_collection.load_from_file()

    for scan in scan_mapping:
        spectrum_obj = spectrum_collection.scandict[int(scan)]
        #Lets determine if the strings are actually ambiguous
        ambiguous_list = ming_ambiguity_library.collapse_ambiguous_from_annotations_list(scan_mapping[scan])
        #print(ambiguous_list)
        if len(ambiguous_list) == 1:
            score_summary = {}
            score_summary["ambiguity_total_score"] = -1
            score_summary["first_unique_count"] = -1
            score_summary["second_unique_count"] = -1
            score_summary["first_unique_intensity"] = -1
            score_summary["second_unique_intensity"] = -1
            score_summary["first_second_unique_ratio"] = -1

            return_ambiguity_mapping[scan] = score_summary

            continue

        if len(ambiguous_list) > 2:
            score_summary = {}
            score_summary["ambiguity_total_score"] = 10
            score_summary["first_unique_count"] = 10
            score_summary["second_unique_count"] = 10
            score_summary["first_unique_intensity"] = 10
            score_summary["second_unique_intensity"] = 10
            score_summary["first_second_unique_ratio"] = -1

            return_ambiguity_mapping[scan] = score_summary
            continue

        peptide_to_extracted_peaks_mapping = {}
        for peptide in ambiguous_list:
            theoreteical_peaks = ming_psm_library.create_theoretical_peak_map(peptide, ["b", "y"])
            original_peaks = spectrum_obj.peaks
            extracted_peaks = extract_annotated_peaks(theoreteical_peaks, original_peaks, peak_tolerance)
            peptide_to_extracted_peaks_mapping[peptide] = extracted_peaks

            #print("Original:\t%d\tExtracted:\t%d" % (len(original_peaks), len(extracted_peaks)))
            #print(original_peaks)
            #print(extracted_peaks)
            #print(theoreteical_peaks)

        #Checkout overlap of stuff
        first_peaks = peptide_to_extracted_peaks_mapping[list(peptide_to_extracted_peaks_mapping.keys())[0]]
        second_peaks = peptide_to_extracted_peaks_mapping[list(peptide_to_extracted_peaks_mapping.keys())[1]]
        total_score, reported_alignments = spectrum_alignment.score_alignment(first_peaks, second_peaks, spectrum_obj.mz, spectrum_obj.mz, peak_tolerance)


        first_total = len(first_peaks)
        second_total = len(second_peaks)
        intersection_total = len(reported_alignments)
        first_unique_count = first_total - intersection_total
        second_unique_count = second_total - intersection_total

        #Calculating the explained intensity in each of these
        peaks_1_normed = spectrum_alignment.sqrt_normalize_spectrum(spectrum_alignment.convert_to_peaks(first_peaks))
        peaks_2_normed = spectrum_alignment.sqrt_normalize_spectrum(spectrum_alignment.convert_to_peaks(second_peaks))



        first_aligned_index = []
        second_aligned_index = []

        for alignment in reported_alignments:
            first_aligned_index.append(alignment.peak1)
            second_aligned_index.append(alignment.peak2)

        #intensity values
        first_unique = []
        second_unique = []

        for i in range(len(peaks_1_normed)):
            if not i in first_aligned_index:
                first_unique.append(peaks_1_normed[i][1])

        for i in range(len(peaks_2_normed)):
            if not i in second_aligned_index:
                second_unique.append(peaks_2_normed[i][1])

        first_unique_intensity = sum(i[0] * i[1] for i in zip(first_unique, first_unique))
        second_unique_intensity = sum(i[0] * i[1] for i in zip(second_unique, second_unique))

        first_second_unique_ratio = 0
        try:
            first_second_unique_ratio = min(first_unique_intensity, second_unique_intensity) / max(first_unique_intensity, second_unique_intensity)
        except KeyboardInterrupt:
            raise
        except:
            first_second_unique_ratio = 10

        if first_second_unique_ratio > 10:
            first_second_unique_ratio = 10

        #print(reported_alignments)
        #print(peaks_1_normed)
        #print("FirstCount\t%d\tSecondCount\t%d\tFirstInt\t%f\tSecondInt\t%f" % (first_unique_count, second_unique_count, first_unique_intensity, second_unique_intensity))

        score_summary = {}
        score_summary["ambiguity_total_score"] = total_score
        score_summary["first_unique_count"] = first_unique_count
        score_summary["second_unique_count"] = second_unique_count
        score_summary["first_unique_intensity"] = first_unique_intensity
        score_summary["second_unique_intensity"] = second_unique_intensity
        score_summary["first_second_unique_ratio"] = first_second_unique_ratio

        return_ambiguity_mapping[scan] = score_summary

    return return_ambiguity_mapping


def main():
    psms_input_file = sys.argv[1]
    input_spectrum_folder = sys.argv[2]
    output_psms_file = sys.argv[3]

    psms_row, psm_table = ming_fileio_library.parse_table_with_headers(psms_input_file)

    peak_tolerance = 0.1

    #Determine which ones have possible bad ambiguity
    spectrum_to_number_psms_dict = defaultdict(lambda: 0)
    psm_table["spectrum_unique_key"] = []
    for i in range(psms_row):
        filename = psm_table["filename"][i]
        scan = psm_table["scan"][i]
        key = filename + ":" + scan
        psm_table["spectrum_unique_key"].append(key)
        spectrum_to_number_psms_dict[key] += 1



    psm_table["numberpsms"] = []
    spectra_to_reconsider = defaultdict(lambda: defaultdict(list))
    for i in range(psms_row):
        filename = psm_table["filename"][i]
        scan = psm_table["scan"][i]
        key = filename + ":" + scan
        number_of_psms_per_spectrum = spectrum_to_number_psms_dict[key]
        psm_table["numberpsms"].append(number_of_psms_per_spectrum)

        if number_of_psms_per_spectrum > 1:
            spectra_to_reconsider[filename][scan].append(psm_table["sequence"][i][:-2])

    spectrum_to_ambiguity_mapping = {}
    for filename in spectra_to_reconsider:
        scan_mapping = spectra_to_reconsider[filename]
        parameter_object = {}
        parameter_object["filename"] = os.path.join(input_spectrum_folder, filename)
        parameter_object["scan_mapping"] = scan_mapping
        scan_ambiguity_mapping = calculated_ambiguity(parameter_object, peak_tolerance)
        for key in scan_ambiguity_mapping:
            full_spectrum_key = "%s:%s" % (filename, key)
            spectrum_to_ambiguity_mapping[full_spectrum_key] = scan_ambiguity_mapping[key]

    psm_table["ambiguity_total_score"] = []
    psm_table["first_unique_count"] = []
    psm_table["second_unique_count"] = []
    psm_table["first_unique_intensity"] = []
    psm_table["second_unique_intensity"] = []
    psm_table["first_second_unique_ratio"] = []
    for i in range(psms_row):
        filename = psm_table["filename"][i]
        scan = psm_table["scan"][i]
        key = filename + ":" + scan
        if key in spectrum_to_ambiguity_mapping:
            psm_table["ambiguity_total_score"].append(spectrum_to_ambiguity_mapping[key]["ambiguity_total_score"])
            psm_table["first_unique_count"].append(spectrum_to_ambiguity_mapping[key]["first_unique_count"])
            psm_table["second_unique_count"].append(spectrum_to_ambiguity_mapping[key]["second_unique_count"])
            psm_table["first_unique_intensity"].append(spectrum_to_ambiguity_mapping[key]["first_unique_intensity"])
            psm_table["second_unique_intensity"].append(spectrum_to_ambiguity_mapping[key]["second_unique_intensity"])
            psm_table["first_second_unique_ratio"].append(spectrum_to_ambiguity_mapping[key]["first_second_unique_ratio"])
        else:
            psm_table["ambiguity_total_score"].append(-1)
            psm_table["first_unique_count"].append(-1)
            psm_table["second_unique_count"].append(-1)
            psm_table["first_unique_intensity"].append(-1)
            psm_table["second_unique_intensity"].append(-1)
            psm_table["first_second_unique_ratio"].append(-1)


    ming_fileio_library.write_dictionary_table_data(psm_table, output_psms_file)




if __name__ == "__main__":
    main()
