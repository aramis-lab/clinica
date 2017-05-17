# -*- coding: utf-8 -*-
"""
    Created on Mon May 15 18:43:59 2017
    
    @author: Pascal
"""

import numpy as np

apoe_snp = ['rs429358', 'rs7412']

# put in BIM file  (position starts with 0)
# 19	rs769451	71.5539	50102751	G	T (563951)
# 19	rs429358	71.5558	50103781	T	C (to be added)
# 19	rs7412	71.5561	50103919	T	C (to be added)
# 19	rs11666145	71.5595	50105714	0	C (563952)

#############################
# Add APOE for BED files

# def import_apoe(file_path):
#     apoe_file = np.genfromtxt(file_path, delimiter=',', skip_header=1, dtype = None, names = ['Phase','ID','RID','SITEID','VISCODE','USERDATE','USERDATE2','APTESTDT','APGEN1','APGEN2','APVOLUME','APRECEIVE','APAMBTEMP','APRESAMP','APUSABLE','update_stamp'])
#     dict = {}
#     for i in range(len(apoe_file)):
#         if apoe_file['Phase'][i] == 'ADNI1':
#             #print apoe_file['APGEN1'][i], apoe_file['APGEN2'][i]
#             dict[apoe_file['RID'][i]] = count_minor_allele_apoe(apoe_file['APGEN1'][i], apoe_file['APGEN2'][i])
#     return dict
#
#
# def reindex_apoe(dict, patient_list, adni_merge):
#     num_patients = len(patient_list)
#     apoe_snp_patient = np.zeros((num_patients, 2))
#     for i in range(num_patients):
#         patient = patient_list[i]
#         rid = int(adni_merge['RID'][np.where(adni_merge['PTID'] == patient)[0][0]])
#         apoe_snp_patient[i] = dict[rid]
#     return apoe_snp_patient


# def count_minor_allele_apoe(epsilon_rs429358, epsilon_rs7412):
#     # rs429358 (T common); rs7412 (T common)
#     if (epsilon_rs429358 == 2) & (epsilon_rs7412 == 2):
#         #Apo2/2 gs268 (T;T) (T;T)
#         count = [0, 0]
#     elif (epsilon_rs429358 == 2) & (epsilon_rs7412 == 3):
#         #Apo2/3 gs269 (T;T) (C;T)
#         count = [0, 1]
#     elif (epsilon_rs429358 == 2) & (epsilon_rs7412 == 4):
#         #Apo2/4 gs270 (C;T) (C;T) ambiguous with Apo1/3
#         count = [1, 1]
#     elif (epsilon_rs429358 == 3) & (epsilon_rs7412 == 3):
#         #Apo3/3 gs246 (T;T) (C;C) the most common
#         count = [0, 2]
#     elif (epsilon_rs429358 == 3) & (epsilon_rs7412 == 4):
#         #Apo3/4 gs141 (C;T) (C;C)
#         count = [1, 2]
#     elif (epsilon_rs429358 == 4) & (epsilon_rs7412 == 4):
#         #Apo4/4 gs216 (C;C) (C;C) ~11x increased Alzheimer's risk
#         count = [2, 2]
#     else:
#         count = []
#         print 'problem with', epsilon_rs429358, epsilon_rs7412
#     return count


#############################
# Add APOE for PED files

def import_apoe_ACGT(file_path):
    apoe_file = np.genfromtxt(file_path, delimiter=',', skip_header=1, dtype = None, names = ['Phase','ID','RID','SITEID','VISCODE','USERDATE','USERDATE2','APTESTDT','APGEN1','APGEN2','APVOLUME','APRECEIVE','APAMBTEMP','APRESAMP','APUSABLE','update_stamp'])
    dict = {}
    for i in range(len(apoe_file)):
        if apoe_file['Phase'][i] == 'ADNI1':
            #print apoe_file['APGEN1'][i], apoe_file['APGEN2'][i]
            dict[apoe_file['RID'][i]] = count_minor_allele_apoe_ACGT(apoe_file['APGEN1'][i], apoe_file['APGEN2'][i])
    return dict


def reindex_apoe_ACGT(dict, patient_list, adni_merge):
    num_patients = len(patient_list)
    apoe_snp_patient = np.zeros((num_patients, 4), dtype = '|S10')
    for i in range(num_patients):
        patient = patient_list[i]
        rid = int(adni_merge['RID'][np.where(adni_merge['PTID'] == patient)[0][0]])
        apoe_snp_patient[i] = dict[rid]
    return apoe_snp_patient


def count_minor_allele_apoe_ACGT(epsilon_rs429358, epsilon_rs7412):
    # rs429358 (T common); rs7412 (T common)
    if (epsilon_rs429358 == 2) & (epsilon_rs7412 == 2):
        #Apo2/2 gs268 (T;T) (T;T)
        sequence = ['T', 'T', 'T', 'T']
    elif (epsilon_rs429358 == 2) & (epsilon_rs7412 == 3):
        #Apo2/3 gs269 (T;T) (C;T)
        sequence = ['T', 'T', 'C', 'T']
    elif (epsilon_rs429358 == 2) & (epsilon_rs7412 == 4):
        #Apo2/4 gs270 (C;T) (C;T) ambiguous with Apo1/3
        sequence = ['C', 'T', 'C', 'T']
    elif (epsilon_rs429358 == 3) & (epsilon_rs7412 == 3):
        #Apo3/3 gs246 (T;T) (C;C) the most common
        sequence = ['T', 'T', 'C', 'C']
    elif (epsilon_rs429358 == 3) & (epsilon_rs7412 == 4):
        #Apo3/4 gs141 (C;T) (C;C)
        sequence = ['C', 'T', 'C', 'C']
    elif (epsilon_rs429358 == 4) & (epsilon_rs7412 == 4):
        #Apo4/4 gs216 (C;C) (C;C) ~11x increased Alzheimer's risk
        sequence = ['C', 'C', 'C', 'C']
    else:
        sequence = []
        print 'problem with', epsilon_rs429358, epsilon_rs7412
    return sequence

