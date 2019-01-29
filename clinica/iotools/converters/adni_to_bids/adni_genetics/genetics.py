# coding: utf8

"""
    Generate new BED/BIM/FAM files by adding APOE and changing patients ID.

    Since APOE has not been genotyped in the original BED/BIM files, but is one
    of the main genes that is responsible for AD, we decided to convert also
    the genetics files of ADNI adding this information during the conversion.

    All the input genetic files .BED, .FAM and .BIM will be stored in a folder
    called Genetics on the top level of
     the BIDS folder, with the following modifications:
     - BED: original file plus APOE
     - FAM: the second column (that contain the original subject id) is
       replaced with the bids id and the fourth and the
      fifth columns ( sex and affection) are replaced with 0
     - BIM: original file plus APOE

    The information regarding APOE will be extracted from the  ADNI clinical
    file APOERS.tsv.
    The output name chosen for the BIDS converted files is
    "ADNI_1_GWAS_Human610-Quad_PlusAPOE".
"""

__author__ = "Pascal Lu"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Sabrina Fontanella"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Pascal Lu"
__email__ = "pascal.lu@icm-institute.org"
__status__ = "Completed"


def convert_genetics(clinical_data_path, input_gen_dir, dest_dir, plink_file_name='ADNI_cluster_01_forward_757LONI',
                     output_file_name='ADNI_1_GWAS_Human610-Quad_PlusAPOE'):
    """
    Convert the genetics data of ADNI to BIDS specification with the
    modification explained above

    Args:
        clinical_data_path: path to the clinical data

        input_gen_dir: path to the folder with the original genetic files
            (.BED, .FAM and .BIM)

        dest_dir: path to the bids directory

        plink_file_name: name of the plink file (if the user don't modify the
            original files downloaded from the website
            is not necessary to modify this variable)

        output_file_name: output name of the converted files

    """
    import os

    path_to_adnimerge = os.path.join(clinical_data_path, 'ADNIMERGE.csv')
    path_to_apoe = os.path.join(clinical_data_path, 'APOERES.csv')
    path_to_plink_files = os.path.join(input_gen_dir, plink_file_name)
    path_to_new_plink_files = os.path.join(dest_dir, 'genetics',
                                           output_file_name)

    if not os.path.exists(os.path.join(dest_dir, 'genetics')):
        os.makedirs(os.path.join(dest_dir, 'genetics'))

    update_plink_files(path_to_adnimerge, path_to_apoe, path_to_plink_files, path_to_new_plink_files)


def update_plink_files(path_to_adnimerge, path_to_apoe, path_to_plink_files, path_to_new_plink_files):
    """
    Add the genotype of APOE and change the patients' ID

    Args:

        path_to_adnimerge: path to the original ADNIMERGE.csv file of ADNI

        path_to_apoe: path the original APOERES.csv file of ADNI

        path_to_plink_files: path to the original folder containing the bed,

        bim and fam file

        path_to_new_plink_files: path to the new plink folder

    """
    import os

    # convert bed bim fam to ped map nof
    os.system('plink --bfile ' + path_to_plink_files + ' --recode --out ' + path_to_plink_files)

    # add apoe to the ped and map files
    add_apoe_ped_files(path_to_adnimerge, path_to_apoe, path_to_plink_files, path_to_new_plink_files)

    # convert ped map nof to bed bim fam
    os.system('plink --file ' + path_to_new_plink_files + ' --make-bed --out ' + path_to_new_plink_files)

    # update fam files by changing patients' ID
    update_ID_fam(path_to_new_plink_files)

    # remove ped, map, nof, hh, log files
    os.remove(path_to_plink_files + '.map')
    if os.path.exists(path_to_plink_files + '.nof'):
        os.remove(path_to_plink_files + '.nof')
    os.remove(path_to_plink_files + '.ped')
    os.remove(path_to_plink_files + '.hh')
    os.remove(path_to_plink_files + '.log')
    os.remove(path_to_new_plink_files + '.map')
    if os.path.exists(path_to_new_plink_files + '.nof'):
        os.remove(path_to_new_plink_files + '.nof')
    os.remove(path_to_new_plink_files + '.ped')
    os.remove(path_to_new_plink_files + '.hh')
    os.remove(path_to_new_plink_files + '.log')


def update_ID_fam(path_to_plink_files):
    """
    Change the patients' ID, replace sex value by 0 and replace affection value by 0

    Args:
        path_to_plink_files: path the plink file

    """
    import csv
    import sys
    import numpy as np

    # load fam file
    csv.field_size_limit(sys.maxsize)
    cr = csv.reader(open(path_to_plink_files + ".fam", "rb"))
    tab = []
    for row in cr:
        tab.append(row[0].split(" "))

    patient_table = np.array(tab)  # its shape is (757, 6)
    new_patient_table = np.zeros(patient_table.shape, dtype='|S20')
    for i in range(len(patient_table)):
        new_patient_table[i, 0] = patient_table[i, 0]
        new_patient_table[i, 1] = 'sub-ADNI' + patient_table[i, 1].replace('_', '')  # change ID
        new_patient_table[i, 2] = patient_table[i, 2]
        new_patient_table[i, 3] = patient_table[i, 3]
        new_patient_table[i, 4] = '0'  # replace sex by 0
        new_patient_table[i, 5] = '0'  # replace affection by 0

    # then write to new fam file
    with open(path_to_plink_files + ".fam", "wb") as stream:
        c = csv.writer(stream)
        for i in range(len(patient_table)):
            l_new_pat = new_patient_table[i, :].tolist()
            c.writerow([' '.join(l_new_pat)])
    print "subjects' ID in fam file has changed."


def add_apoe_ped_files(path_to_adnimerge, path_to_apoe, path_to_plink_files, path_to_new_plink_files):
    """
   Add the genotype of APOE inside the DNA sequence

    Args:

        path_to_adnimerge: path to the original ADNIMERGE.csv file of ADNI

        path_to_apoe: path the original APOERES.csv file of ADNI

        path_to_plink_files: path to the old plink file to modify

        path_to_new_plink_files: path to the new plink file modified

    Returns:

    """
    import csv
    import sys
    import numpy as np
    import pandas as pd

    # load ped file
    csv.field_size_limit(sys.maxsize)
    cr = csv.reader(open(path_to_plink_files + ".ped", "rb"))
    tab = []
    for row in cr:
        tab.append(row[0].split(" "))
    snp = np.array(tab)  # its shape is (757, 1241808)
    patient_list = snp[:, 1]

    # add apoe to ped file
    adni_merge = pd.io.parsers.read_csv(path_to_adnimerge, sep=',')
    dict_patients_apoe = import_apoe_ACGT(path_to_apoe)
    apoe_snp_patient = reindex_apoe_ACGT(dict_patients_apoe, patient_list, adni_merge)
    final_snp = np.concatenate((snp[:, :(563952 * 2 + 6)], apoe_snp_patient, snp[:, (563952 * 2 + 6):]), axis=1)

    # write to ped file
    print 'ped file shape: ', final_snp.shape
    with open(path_to_new_plink_files + ".ped", "wb") as stream:
        c = csv.writer(stream)
        for r in range(len(final_snp)):
            l_new_pat = final_snp[r, :].tolist()
            c.writerow([' '.join(l_new_pat)])
    print 'new ped file created'

    # load map file and write new map file
    cr2 = csv.reader(open(path_to_plink_files + ".map", "rb"))
    tab2 = []
    for row2 in cr2:
        tab2.append(row2[0].split("\t"))
    snp_list = np.array(tab2)

    # put in BIM file  (position starts with 0)
    # 19	rs769451	71.5539	50102751	G	T (563951)
    # 19	rs429358	71.5558	50103781	T	C (to be added)
    # 19	rs7412	71.5561	50103919	T	C (to be added)
    # 19	rs11666145	71.5595	50105714	0	C (563952)

    final_snp_list = np.concatenate((snp_list[:563952],
                                     [['19', 'rs429358', '71.5558', '50103781']],
                                     [['19', 'rs7412', '71.5561', '50103919']],
                                     snp_list[563952:]), axis=0)

    with open(path_to_new_plink_files + ".map", "wb") as stream:
        c = csv.writer(stream)
        for r in range(len(final_snp_list)):
            l_new_pat = final_snp_list[r, :].tolist()
            c.writerow(['\t'.join(l_new_pat)])
    print 'new map file created'


def import_apoe_ACGT(file_path, bed=0):
    """
    For all subjects, it returns the genotype of APOE

    Args:

        file_path: path to the apoe file

        bed: 0 or 1. If 0, the genotype will be seen as a DNA sequence. If 1,
            the genotype will be a sequence of numbers (0,1,2), where the number of
            minor variants is indicated for each snp

    Returns:

        dict: dictonary where the keys are the subjects id and the values are
            the sequence

    """
    import pandas as pd
    apoe_file = pd.io.parsers.read_csv(file_path, sep=',')
    dict = {}
    for i in range(len(apoe_file)):
        if apoe_file['Phase'][i] == 'ADNI1':
            seq, count = count_minor_allele_apoe_ACGT(apoe_file['APGEN1'][i], apoe_file['APGEN2'][i])
            if bed:
                dict[apoe_file['RID'][i]] = count
            else:
                dict[apoe_file['RID'][i]] = seq
    return dict


def reindex_apoe_ACGT(dict, patient_list, adni_merge, bed=0):
    """
    Create a table with the DNA sequence of gene APOE ordered by the list of patients

    Args:
        dict: the output of the method import_apoe_ACGT

        patient_list: list of patients

        adni_merge: path to the ADNIMERGE.csv file of ADNI

        bed: 0 or 1. If 0, the genotype will be seen as a DNA sequence. If 1,
            the genotype will be a sequence of numbers (0,1,2), where the
            number of minor variants is indicated for each snp. It must be the
            same value passed to the method import_apoe_ACGT

    Returns:
        apoe_snp_patient: a table of DNA sequence for the gene APOE ordered by
            the list of patients
    """
    import numpy as np
    num_patients = len(patient_list)
    if bed:
        apoe_snp_patient = np.zeros((num_patients, 2))
    else:
        apoe_snp_patient = np.zeros((num_patients, 4), dtype='|S10')
    for i in range(num_patients):
        patient = patient_list[i]
        rid = int(adni_merge['RID'][np.where(adni_merge['PTID'] == patient)[0][0]])
        apoe_snp_patient[i] = dict[rid]
    return apoe_snp_patient


def count_minor_allele_apoe_ACGT(epsilon_first_parent, epsilon_second_parent):
    """
    Convert the allele to the DNA sequence

    If the APOE allele is epsilon4/epsilon4, this method will return the
    sequence C,C,C,C

    Args:

        epsilon_first_parent: allele from the first parent

        epsilon_second_parent: allele from the second parent

    Returns:

        sequence: the two first letters correspond to snp rs429358, the two
            last letters correspond to snp rs7412

        count: the first number correspond to the number of minor variants for
            rs429358, the second number correspond to the number of minor
            variants for rs7412
    """
    # rs429358 (T common); rs7412 (T common)
    if (epsilon_first_parent == 2) & (epsilon_second_parent == 2):
        # Apo2/2 gs268 (T;T) (T;T)
        sequence = ['T', 'T', 'T', 'T']
        count = [0, 0]
    elif (epsilon_first_parent == 2) & (epsilon_second_parent == 3):
        # Apo2/3 gs269 (T;T) (C;T)
        sequence = ['T', 'T', 'C', 'T']
        count = [0, 1]
    elif (epsilon_first_parent == 2) & (epsilon_second_parent == 4):
        # Apo2/4 gs270 (C;T) (C;T) ambiguous with Apo1/3
        sequence = ['C', 'T', 'C', 'T']
        count = [1, 1]
    elif (epsilon_first_parent == 3) & (epsilon_second_parent == 3):
        # Apo3/3 gs246 (T;T) (C;C) the most common
        sequence = ['T', 'T', 'C', 'C']
        count = [0, 2]
    elif (epsilon_first_parent == 3) & (epsilon_second_parent == 4):
        # Apo3/4 gs141 (C;T) (C;C)
        sequence = ['C', 'T', 'C', 'C']
        count = [1, 2]
    elif (epsilon_first_parent == 4) & (epsilon_second_parent == 4):
        # Apo4/4 gs216 (C;C) (C;C) ~11x increased Alzheimer's risk
        sequence = ['C', 'C', 'C', 'C']
        count = [2, 2]
    else:
        sequence = []
        count = []
        print 'missing combination for ', epsilon_first_parent, epsilon_second_parent
    return sequence, count
