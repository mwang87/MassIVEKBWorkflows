#!/usr/bin/python

import re
import distance

#Return 0 for not ambiguous, 1 for ambiguous
def collapse_ambiguous_from_annotations_list(list_of_peptides):
    new_peptide_list = []
    for peptide in list_of_peptides:
        new_peptide = peptide.replace("I", "L")
        new_peptide_list.append(new_peptide)

    #Let make sure its not just I/L substitutions
    new_peptide_list = list(set(new_peptide_list))
    return new_peptide_list

def categorize_peptide_distance(annotation1, annotation2):
    #Determining if it is I/L Substitution
    if annotation1.replace("I", "L") == annotation2.replace("I", "L"):
        #I/L Substitution
        return "I/L Substitution"

    annotation1_sequence_only = re.sub(r'[0-9.+-]+', '', annotation1)
    annotation2_sequence_only = re.sub(r'[0-9.+-]+', '', annotation2)

    string_distance = distance.nlevenshtein(annotation1_sequence_only, annotation2_sequence_only, method=1)
    #Detecting Site locatization of PTMs
    if string_distance < 0.01:
        return "PTM Localization"

    hamming_distance = 0

    if len(annotation1_sequence_only) == len(annotation2_sequence_only):
        hamming_distance = distance.hamming(annotation1_sequence_only, annotation2_sequence_only)

        if hamming_distance == 2:
            return "Double Amino Substitution"

        if hamming_distance == 1:
            #Seeing if it is a deamidation
            annotation1_contains_deamidation = False
            annotation2_contains_deamidation = False

            if annotation1.find("+0.984") != -1:
                annotation1_contains_deamidation = True
            if annotation2.find("+0.984") != -1:
                annotation2_contains_deamidation = True

            if annotation1_contains_deamidation != annotation2_contains_deamidation:
                #Probably should check for Q->E
                return "Deamidation"


            #Checking for Q->K Substitution

    #Determining String Distance
    string_distance = distance.nlevenshtein(annotation1, annotation2, method=1)

    return "UNKNOWN"
