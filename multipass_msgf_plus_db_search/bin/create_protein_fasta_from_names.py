#!/usr/bin/python


import sys
import getopt
import os
import ming_protein_library
import ming_psm_library
import ming_parallel_library
import library_creation
import json

def usage():
    print "<input fasta> <input protein names> <output fasta file>"



def main():
    input_fasta_filename = sys.argv[1]
    input_protein_names_file = sys.argv[2]
    output_fasta_file = sys.argv[3]

    proteome = ming_protein_library.parse_fasta_proteome_file(input_fasta_filename)

    protein_list = json.loads(open(input_protein_names_file, "r").read())
    protein_list = list(set(protein_list))

    new_protein_list = []
    decoy_prot_count = 0
    for protein in protein_list:
        if protein.find("REV_") != -1:
            new_protein_list.append(protein[4:])
            decoy_prot_count += 1
        else:
            new_protein_list.append(protein)
    print "Decoy Prot Count: " + str(decoy_prot_count)
    print "Old prot list len: "  + str(len(protein_list))
    print "New prot list len: "  + str(len(new_protein_list))
    new_protein_list = list(set(new_protein_list))
    print "New prot list len: "  + str(len(new_protein_list))

    output_fasta_string = ""
    for protein in new_protein_list:
        output_fasta_string += proteome.get_protein(protein).to_fasta()

    open(output_fasta_file, "w").write(output_fasta_string)

if __name__ == "__main__":
    main()
