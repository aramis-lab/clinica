# -*- coding: utf-8 -*-
"""
    Created on Mon May 15 18:43:59 2017
    
    @author: Pascal
"""

import numpy as np

# def import_genetics(path_to_plink_files):
#     from plinkio import plinkfile
#     """
#     Input:
#         - path to a folder containing a bed, bim and fam file
#     -------------------------------------------------------
#     Output:
#         - snp_list : list of snps
#         - patients: list of patients
#         - matrix with patient as line and snp as column
#     -------------------------------------------------------
#     Note: how to read the bed/bim/fam files:
#         ## affection
#         - unknown: 0 (fam) -> -9 (python)
#         - unaffected: 1 (fam) -> 0 (python)
#         - affected: 2 (fam) -> 1 (python)
#
#         ## sex
#         - male: 1 (fam) -> 0 (python)
#         - female: 2 (fam) -> 1 (python)
#
#         ## genotype
#         - genotype 0: code 00 Homozygote "0"/"0"
#         - genotype 1: code 01 Heterozygote
#         - genotype 2: code 11 Homozygote "1"/"1"
#         - genotype 3: unknown
#     -------------------------------------------------------
#     More details are available on https://web.njit.edu/~zhiwei/GAS_101.pdf
#     and on https://github.com/mfranberg/libplinkio
#     and on http://www.gwaspi.org/?page_id=671
#     """
#
#     plink_file = plinkfile.open(path_to_plink_files)
#     if not plink_file.one_locus_per_row():
#         print("This script requires that snps are rows and samples columns.")
#         exit(1)
#
#     sample_list = plink_file.get_samples()
#     locus_list = plink_file.get_loci()
#
#     """
#     for sample in sample_list:
#         print sample.fid, sample.iid, sample.father_iid, sample.mother_iid, sample.sex, sample.affection
#     for locus in locus_list:
#         print locus.chromosome, locus.name, locus.position, locus.bp_position, locus.allele1, locus.allele2
#     """
#
#     num_snp = len(locus_list)
#     num_patient = len(sample_list)
#
#     matrix_patient_snp = np.zeros((num_snp, num_patient))
#     patient_list = np.zeros(num_patient, dtype='|S30')
#     snp_list = np.zeros(num_snp, dtype='|S30')
#     for i, row, locus in zip( range(num_snp), plink_file, locus_list):
#         matrix_patient_snp[i] = row
#         snp_list[i] = locus.name
#     matrix_patient_snp = matrix_patient_snp.T
#     for i, sample in zip(range(num_patient), sample_list):
#         patient_list[i] = sample.iid
#     return snp_list, patient_list, matrix_patient_snp



def import_adnimerge(path_to_adnimerge):
    """
    Input:
        - path to ADNIMERGE.csv
    -------------------------------------------------------
    Output:
        - data in ADNIMERGE.csv
    """
    adni_file = np.genfromtxt(path_to_adnimerge, delimiter=';', skip_header=0, dtype = None, names = ['RID', 'PTID', 'VISCODE', 'SITE', 'COLPROT', 'ORIGPROT', 'EXAMDATE', 'DX_bl', 'AGE', 'PTGENDER', 'PTEDUCAT', 'PTETHCAT', 'PTRACCAT', 'PTMARRY', 'APOE4', 'FDG', 'PIB', 'AV45', 'CDRSB', 'ADAS11', 'ADAS13', 'MMSE', 'RAVLT_immediate', 'RAVLT_learning', 'RAVLT_forgetting', 'RAVLT_perc_forgetting', 'FAQ', 'MOCA', 'EcogPtMem', 'EcogPtLang', 'EcogPtVisspat', 'EcogPtPlan', 'EcogPtOrgan', 'EcogPtDivatt', 'EcogPtTotal', 'EcogSPMem', 'EcogSPLang', 'EcogSPVisspat', 'EcogSPPlan', 'EcogSPOrgan', 'EcogSPDivatt', 'EcogSPTotal', 'Ventricles', 'Hippocampus', 'WholeBrain', 'Entorhinal', 'Fusiform', 'MidTemp', 'ICV', 'DX', 'EXAMDATE_bl', 'CDRSB_bl', 'ADAS11_bl', 'ADAS13_bl', 'MMSE_bl', 'RAVLT_immediate_bl', 'RAVLT_learning_bl', 'RAVLT_forgetting_bl', 'RAVLT_perc_forgetting_bl', 'FAQ_bl', 'Ventricles_bl', 'Hippocampus_bl', 'WholeBrain_bl', 'Entorhinal_bl', 'Fusiform_bl', 'MidTemp_bl', 'ICV_bl', 'MOCA_bl', 'EcogPtMem_bl', 'EcogPtLang_bl', 'EcogPtVisspat_bl', 'EcogPtPlan_bl', 'EcogPtOrgan_bl', 'EcogPtDivatt_bl', 'EcogPtTotal_bl', 'EcogSPMem_bl', 'EcogSPLang_bl', 'EcogSPVisspat_bl', 'EcogSPPlan_bl', 'EcogSPOrgan_bl', 'EcogSPDivatt_bl', 'EcogSPTotal_bl', 'FDG_bl', 'PIB_bl', 'AV45_bl', 'Years_bl', 'Month_bl', 'Month', 'M', 'update_stamp'])
    return adni_file



def import_preselected_adnimerge(path_to_adnimerge, patient_list):
    """
    Input:
        - path to ADNIMERGE.csv
        - patient_list
    -------------------------------------------------------
    Output:
        - data in ADNIMERGE.csv with patients that are in patient_list
    """
    adni_file = import_adnimerge(path_to_adnimerge)
    preselected_adni_file = np.zeros((adni_file.shape), dtype = (adni_file).dtype)
    num_lines = 0
    
    for i in range(len(adni_file)):
        if (adni_file[i]['PTID'] in patient_list):
            preselected_adni_file[num_lines] = adni_file[i]
            num_lines += 1

    return preselected_adni_file[:num_lines]
