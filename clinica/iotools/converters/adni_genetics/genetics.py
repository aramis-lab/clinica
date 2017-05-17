def convert_genetics(clinical_data_path, genetics_files_paths, dest_dir, plink_file_name = 'ADNI_cluster_01_forward_757LONI'):
    import os

    output_file_name = 'ADNI_1_GWAS_Human610-Quad_PlusAPOE'



    path_to_adnimerge = os.path.join(clinical_data_path, 'ADNIMERGE.csv')
    path_to_apoe = os.path.join(clinical_data_path, 'APOERES.csv')
    path_to_plink_files = os.path.join(genetics_files_paths, plink_file_name)
    path_to_new_plink_files = os.path.join(dest_dir, 'genetics',
                                           output_file_name)

    if not os.path.exists(os.path.join(dest_dir, 'genetics')):
        os.makedirs(os.path.join(dest_dir, 'genetics'))

    update_plink_files(path_to_adnimerge, path_to_apoe, path_to_plink_files, path_to_new_plink_files)


def update_plink_files(path_to_adnimerge, path_to_apoe, path_to_plink_files, path_to_new_plink_files):
    import os
    import shutil

    # convert bed bim fam to ped map nof
    os.system('plink --bfile ' + path_to_plink_files + ' --recode --out ' + path_to_plink_files)

    # add apoe to the ped and map files
    add_apoe_ped_files(path_to_adnimerge, path_to_apoe, path_to_plink_files, path_to_new_plink_files)

    # copy nof file without modification
    shutil.copyfile(path_to_plink_files + '.nof', path_to_new_plink_files + '.nof')

    # convert ped map nof to bed bim fam
    os.system('plink --file ' + path_to_new_plink_files + ' --make-bed --out ' + path_to_new_plink_files)

    # update fam files by changing patients' ID
    update_ID_fam(path_to_new_plink_files)

    # remove ped, map, nof, hh, log files
    os.remove(path_to_plink_files + '.map')
    os.remove(path_to_plink_files + '.nof')
    os.remove(path_to_plink_files + '.ped')
    os.remove(path_to_plink_files + '.hh')
    os.remove(path_to_plink_files + '.log')
    os.remove(path_to_new_plink_files + '.map')
    os.remove(path_to_new_plink_files + '.nof')
    os.remove(path_to_new_plink_files + '.ped')
    os.remove(path_to_new_plink_files + '.hh')
    os.remove(path_to_new_plink_files + '.log')


def add_apoe_ped_files(path_to_adnimerge, path_to_apoe, path_to_plink_files, path_to_new_plink_files):
    import csv
    import sys
    import load_data
    import preprocess_apoe
    import numpy as np

    # load ped file
    csv.field_size_limit(sys.maxsize)
    cr = csv.reader(open(path_to_plink_files + ".ped", "rb"))
    tab = []
    for row in cr:
        tab.append(row[0].split(" "))
    snp = np.array(tab)  # its shape is (757, 1241808)
    patient_list = snp[:, 1]

    # add apoe to ped file
    adni_merge = load_data.import_preselected_adnimerge(path_to_adnimerge, patient_list)
    dict_patients_apoe = preprocess_apoe.import_apoe_ACGT(path_to_apoe)
    apoe_snp_patient = preprocess_apoe.reindex_apoe_ACGT(dict_patients_apoe, patient_list, adni_merge)
    final_snp = np.concatenate((snp[:, :(563952 * 2 + 6)], apoe_snp_patient, snp[:, (563952 * 2 + 6):]), axis=1)

    # write to ped file
    print 'ped file shape: ', final_snp.shape
    with open(path_to_new_plink_files + ".ped", "wb") as stream:
        c = csv.writer(stream)
        for r in range(len(final_snp)):
            l = final_snp[r, :].tolist()
            c.writerow([' '.join(l)])
    print 'new ped file created'

    # load map file and write new map file
    cr2 = csv.reader(open(path_to_plink_files + ".map", "rb"))
    tab2 = []
    for row2 in cr2:
        tab2.append(row2[0].split("\t"))
    snp_list = np.array(tab2)
    final_snp_list = np.concatenate((snp_list[:563952],
                                     [['19', 'rs429358', '71.5558', '50103781']],
                                     [['19', 'rs7412', '71.5561', '50103919']],
                                     snp_list[563952:]), axis=0)

    with open(path_to_new_plink_files + ".map", "wb") as stream:
        c = csv.writer(stream)
        for r in range(len(final_snp_list)):
            l = final_snp_list[r, :].tolist()
            c.writerow(['\t'.join(l)])
    print 'new map file created'


def update_ID_fam(path_to_plink_files):
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
        new_patient_table[i, 1] = 'sub_ADNI' + patient_table[i, 1].replace('_', '')  # change ID
        new_patient_table[i, 2] = patient_table[i, 2]
        new_patient_table[i, 3] = patient_table[i, 3]
        new_patient_table[i, 4] = '0'  # replace sex by 0
        new_patient_table[i, 5] = '0'  # replace affection by 0

    # then write to new fam file
    with open(path_to_plink_files + ".fam", "wb") as stream:
        c = csv.writer(stream)
        for i in range(len(patient_table)):
            l = new_patient_table[i, :].tolist()
            c.writerow([' '.join(l)])

    print "subjects' ID in fam file has changed."