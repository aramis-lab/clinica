# coding: utf-8

"""
Load bed/bim/fam files.
"""

__author__ = "Pascal Lu"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Sabrina Fontanella"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Pascal Lu"
__email__ = "pascal.lu@icm-institute.org"
__status__ = "Completed"


def import_genetics(path_to_plink_files):
    """
    Import genetics data

    Note about how to read the bed/bim/fam files:
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

    Args:
        path_to_plink_files: path to a folder containing a bed, bim and fam file

    Returns:
        snp_list: list of snps
        patients: list of patients
        matrix: matrix with patient as line and snp as column

    """
    from plinkio import plinkfile
    import numpy as np

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
    for i, row, locus in zip(range(num_snp), plink_file, locus_list):
        matrix_patient_snp[i] = row
        snp_list[i] = locus.name
    matrix_patient_snp = matrix_patient_snp.T
    for i, sample in zip(range(num_patient), sample_list):
        patient_list[i] = sample.iid
    return snp_list, patient_list, matrix_patient_snp
