#!/usr/bin/python

import ming_psm_library
import ming_protein_library

#Checks to see whether this variant is actually consistent within itself
def check_variant_psm_consistency(peptide_variant):
    #Maybe run to a database and grab the spectra
    for psm in peptide_variant.psms:
        print(psm)
    #print "CHECKING THIS SHIT"


def protein_FDR_calculated_protein_FDR(protein_list):
    target_count = 0
    decoy_count = 0

    for protein in protein_list:
        if protein.number_precursors == 0:
            protein.score = -1000

    #sorting descrending score
    protein_list = sorted(protein_list, key=lambda protein: protein.score, reverse=True)

    for protein in protein_list:
        if protein.decoy == 1:
            decoy_count += 1
        else:
            target_count += 1
        fdr = float(decoy_count)/float(decoy_count + target_count)
        protein.fdr = fdr
        if protein.score == -1000:
            protein.fdr = 1.0

    number_of_one_percent_proteins = 0
    for protein in protein_list:
        if protein.fdr < 0.01:
            number_of_one_percent_proteins += 1

    print("Number of proteins under 1%", number_of_one_percent_proteins)

    return protein_list

"""Protein FDR by the method proposed by proteomicsDB. Eliminates the decoy of the target for scoring"""
def protein_FDR_calculated_protein_FDR_protein_elimination(protein_list):
    target_count = 0
    decoy_count = 0

    for protein in protein_list:
        if protein.number_precursors == 0:
            protein.score = -1000

    #sorting descrending score
    protein_list = sorted(protein_list, key=lambda protein: protein.score, reverse=True)

    seen_proteins = set()
    for protein in protein_list:
        if protein.decoy == 1:
            #print("DECOY", protein.protein)
            non_decoy_protein_name = protein.protein.replace("XXX_", "")
            if non_decoy_protein_name in seen_proteins:
                decoy_count = decoy_count
            else:
                decoy_count += 1
        else:
            target_count += 1
        fdr = float(decoy_count)/float(decoy_count + target_count)
        protein.fdr = fdr
        if protein.score == -1000:
            protein.fdr = 1.0

        seen_proteins.add(protein.protein)

    number_of_one_percent_proteins = 0
    for protein in protein_list:
        if protein.fdr < 0.01:
            number_of_one_percent_proteins += 1

    print("Number of proteins under 1%", number_of_one_percent_proteins)

    return protein_list

def protein_FDR_assign_scoring_by_variant_count(proteome, psms, peptide_to_protein_mapping):
    #Setting initial score for proteins
    for protein in proteome.protein_list:
        protein.score = 0
        protein.number_precursors = 0

    for psm in psms:
        annotation = psm.annotation
        score = psm.score
        stripped_sequence = ming_psm_library.strip_sequence(annotation).replace("I", "L")
        protein_list = peptide_to_protein_mapping[stripped_sequence]


        if len(protein_list) > 1:
            continue
            gene_list = []
            for protein in protein_list:
                protein_obj = proteome.protein_map[protein]
                gene_name = protein_obj.gene_name
                gene_list.append(gene_name)
            gene_list = list(set(gene_list))
            if len(gene_list) > 1:
                print(stripped_sequence, len(gene_list), "Skipping because not unique")
                continue


        for protein_name in protein_list:
            if len(protein_name) < 2:
                continue
            proteome.get_protein(protein_name).score += score
            #proteome.get_protein(protein_name).score += 1
            proteome.get_protein(protein_name).number_precursors += 1
