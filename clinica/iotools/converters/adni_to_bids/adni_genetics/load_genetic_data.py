# -*- coding: utf-8 -*-
"""
    Created on Mon May 15 18:43:59 2017
    
    load_data.py load bed/bim/fam files
    
    * It requires plinkio library (https://pypi.python.org/pypi/plinkio)
    * Type for installation: sudo pip install plinkio
    
    @author: Pascal
"""

import numpy as np


def import_genetics(path_to_plink_files):
    from plinkio import plinkfile
    """
    Input:
        - path to a folder containing a bed, bim and fam file
    -------------------------------------------------------
    Output:
        - snp_list : list of snps
        - patients: list of patients
        - matrix with patient as line and snp as column
    -------------------------------------------------------
    Note: how to read the bed/bim/fam files:
        ## affection
        - unknown: 0 (fam) -> -9 (python)
        - unaffected: 1 (fam) -> 0 (python)
        - affected: 2 (fam) -> 1 (python)
        
        ## sex
        - male: 1 (fam) -> 0 (python)
        - female: 2 (fam) -> 1 (python)
        
        ## genotype
        - genotype 0: code 00 Homozygote "0"/"0"
        - genotype 1: code 01 Heterozygote
        - genotype 2: code 11 Homozygote "1"/"1"
        - genotype 3: unknown
    -------------------------------------------------------
    More details are available on https://web.njit.edu/~zhiwei/GAS_101.pdf
    and on https://github.com/mfranberg/libplinkio
    and on http://www.gwaspi.org/?page_id=671
    """
    
    plink_file = plinkfile.open(path_to_plink_files)
    if not plink_file.one_locus_per_row():
        print("This script requires that snps are rows and samples columns.")
        exit(1)
    
    sample_list = plink_file.get_samples()
    locus_list = plink_file.get_loci()

    """
    for sample in sample_list:
        print sample.fid, sample.iid, sample.father_iid, sample.mother_iid, sample.sex, sample.affection
    for locus in locus_list:
        print locus.chromosome, locus.name, locus.position, locus.bp_position, locus.allele1, locus.allele2
    """
    
    num_snp = len(locus_list)
    num_patient = len(sample_list)
    
    matrix_patient_snp = np.zeros((num_snp, num_patient))
    patient_list = np.zeros(num_patient, dtype='|S30')
    snp_list = np.zeros(num_snp, dtype='|S30')
    for i, row, locus in zip( range(num_snp), plink_file, locus_list):
        matrix_patient_snp[i] = row
        snp_list[i] = locus.name
    matrix_patient_snp = matrix_patient_snp.T
    for i, sample in zip(range(num_patient), sample_list):
        patient_list[i] = sample.iid
    return snp_list, patient_list, matrix_patient_snp


if __name__ == '__main__':
    import pandas as pd
    import genetics
    
    path_to_plink_files = '/Users/sabrina.fontanella/ownCloud/genotype/data/ADNI_1_GWAS_Human610-Quad_PlusAPOE/ADNI_1_GWAS_Human610-Quad_PlusAPOE'
    path_to_adnimerge = '/Users/sabrina.fontanella/ownCloud/genotype/code_sabrina/data/ADNIMERGE.csv'
    path_to_apoe = '/Users/sabrina.fontanella/ownCloud/genotype/code_sabrina/data/APOERES.csv'
    path_to_new_plink_files = '/Users/sabrina.fontanella/ownCloud/genotype/code_sabrina/data/genetics/ADNI_1_GWAS_Human610-Quad_PlusAPOE'
    
    
    # load original genetic data
    snp_list, patient_list, matrix_patient_snp = import_genetics(path_to_plink_files)
    print patient_list
    # adni_merge = pd.io.parsers.read_csv(path_to_adnimerge, sep=',')
    # dict_patients_apoe = genetics.import_apoe_ACGT(path_to_apoe, bed=1)
    # apoe_snp_patient = genetics.reindex_apoe_ACGT(dict_patients_apoe, patient_list, adni_merge, bed=1)
    # matrix_patient_snp = np.concatenate((matrix_patient_snp[:, :563952], apoe_snp_patient, matrix_patient_snp[:, 563952:]),axis=1)
    # snp_list = np.concatenate((snp_list[:563952],np.array(['rs429358', 'rs7412']),snp_list[563952:]), axis=0)
    #
    #
    # # load modified genetic data
    # snp_list2, patient_list2, matrix_patient_snp2 = import_genetics(path_to_new_plink_files)
    #
    #
    # # compare and check that the insertion of SNPs belonging to APOE is correct
    # for i in [563951,563952,563953,563954]:
    #     print ''
    #     print 'position : ', i
    #     print 'snp (original file): ', snp_list[i]
    #     print 'snp (new modified file): ', snp_list2[i]
    #     print 'identical data: ', ((matrix_patient_snp[:,i] == matrix_patient_snp2[:,i]).all() | (matrix_patient_snp[:,i] == 2-matrix_patient_snp2[:,i]).all())
