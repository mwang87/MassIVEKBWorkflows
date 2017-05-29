#!/usr/bin/python

import ming_psm_library


#Checks to see whether this variant is actually consistent within itself
def check_variant_psm_consistency(peptide_variant):
    #Maybe run to a database and grab the spectra
    for psm in peptide_variant.psms:
        print(psm)
    #print "CHECKING THIS SHIT"
