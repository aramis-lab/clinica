"""
Covert the HMTC dataset into the BIDS specification.

@author: Jeremy Guillon
"""

from converter_utils import MissingModsTracker, print_statistics
import pkg_resources as pkg
import bids_utils as bids
from shutil import copy
from glob import glob
from os import path
import pandas as pd
import logging
import json
import os


def convert(source_dir, dest_dir):
    """
    Convert the HMTC dataset into the BIDS standard.

    Args:
        source_dir: directory of the input dataset
        dest_dir: output directory
        param: parameters for the conversion
    """

    print "*******************************"
    print "HMTC to BIDS converter"
    print "*******************************"

    # Prepare output directory and log files
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)
    summary_file = open(path.join(dest_dir, 'conversion_summary.txt'), 'w')
    logging.basicConfig(filename=path.join(dest_dir, 'conversion_modalities.log'),
                        format='%(asctime)s %(levelname)s:%(message)s',
                        level=logging.DEBUG,
                        datefmt='%m/%d/%Y %I:%M')

    # Init subjects' information tables
    subject_hmtc_paths = []
    subject_bids_paths = []
    hmtc_ids = []
    bids_ids = []

    # Create the lists of bids_ids extracting the list of the subjects from
    # the dataset
    subject_hmtc_paths = glob(path.join(source_dir, '*'))
    for subject_hmtc_path in subject_hmtc_paths:
        subject_hmtc_id = subject_hmtc_path.split(os.sep)[-1]
        hmtc_ids.append(subject_hmtc_id)
        bids_ids.append('sub-HMTC' + subject_hmtc_id)
        subject_bids_paths.append(path.join(dest_dir, bids_ids[-1]))
        if not os.path.exists(subject_bids_paths[-1]):
            os.makedirs(subject_bids_paths[-1])

    print "Number of subjects: ", len(hmtc_ids)

    df = pd.read_excel(pkg.resource_filename(
        'clinica', 'bids/data/modality_equiv.xlsx'))

    # For each subject extract the list of files and convert them into BIDS
    # specification
    for i in range(len(hmtc_ids)):

        print "Converting: ", subject_hmtc_paths[i]
        logging.info("Converting: " + subject_hmtc_paths[i])
        subject_bids_prefix = bids_ids[i] + "_ses-M00"

        # Create session directory
        session_bids_path = path.join(subject_bids_paths[i], 'ses-M00')
        if not os.path.exists(session_bids_path):
            os.makedirs(session_bids_path)

        # For each expected modality
        for m in range(len(df)):

            mod_hmtc_paths = glob(
                path.join(subject_hmtc_paths[i], df.regex[m]))
            if len(mod_hmtc_paths) > df.expected_num[m]:
                if df.type[m] != 'magnitude' and df.type[m] != 'phasediff':
                    logging.warning('Subject ' + hmtc_ids[i] + ' has ' + str(len(mod_hmtc_paths)) + ' ' +
                                    df.type[m] + ' (' + str(df.expected_num[m]) + ' expected). Only the ' + str(df.expected_num[m]) + ' first one(s) will be copied.')
                    print 'Subject ' + hmtc_ids[i] + ' has ' + str(len(mod_hmtc_paths)) + ' ' + df.type[m] + ' (' + str(df.expected_num[m]) + ' expected). Only the ' + str(df.expected_num[m]) + ' first one(s) will be copied.'
                add_modalities(session_bids_path, subject_bids_prefix, mod_hmtc_paths,
                               df.type[m], df.acq_label[m], df.expected_num[m])
            elif len(mod_hmtc_paths) < df.expected_num[m]:
                logging.warning('Subject ' + hmtc_ids[i] + ' is missing ' +
                                str(df.expected_num[m] - len(mod_hmtc_paths)) + ' ' + df.type[m] + ' file(s). Skipping.')
                print 'Subject ' + hmtc_ids[i] + ' is missing ' + str(df.expected_num[m] - len(mod_hmtc_paths)) + ' ' + df.type[m] + ' file(s). Skipping.'
            else:
                add_modalities(session_bids_path, subject_bids_prefix, mod_hmtc_paths,
                               df.type[m], df.acq_label[m], df.expected_num[m])

        print '---'


def add_modalities(session_dir, filename_prefix, modality_dirs, type, acq_label=None, expected_num=1):

    category_dir = get_category_dir(type)
    if not os.path.exists(path.join(session_dir, category_dir)):
        os.makedirs(path.join(session_dir, category_dir))

    files_ext = get_files_ext(type)
    for file_ext in files_ext:

        n_run = 0

        for modality_dir in modality_dirs:

            n_file = 0
            n_run += 1

            # Check if we have enough files
            if n_run > expected_num and category_dir != 'fmap':
                continue

            files = glob(path.join(modality_dir, '*.' + file_ext))

            if len(files) == 0:
                logging.warning('')  # TODO
            elif len(files) > 1 and category_dir != 'fmap':
                logging.warning('')  # TODO
            else:
                for file in files:
                    n_file += 1
                    filename = filename_prefix
                    if not pd.isnull(acq_label):
                        filename += '_acq-' + str(acq_label)
                    if expected_num > 1 and category_dir != 'fmap':
                        filename += '_run-' + str(n_run)
                    # Special cases of fmaps
                    if type == 'magnitude' and n_run == 1 and expected_num > 1:
                        filename += '_' + type + str(n_file)
                    elif type == 'magnitude' and n_run > 1:
                        continue
                    elif type == 'phasediff' and n_run == 2 and expected_num > 1:
                        filename += '_' + type + str(n_file)
                    elif type == 'phasediff' and n_run < 2:
                        continue
                    else:
                        filename += '_' + type
                    filename += '.' + file_ext
                    copy(file, path.join(session_dir, category_dir, filename))


def get_category_dir(type):
    """Give the BIDS modality folder in function of the modality type

    :param type: the modality type amoung the following list: 'T1w', 'T2w', 'FLAIR', 'bold', 'dwi', 'magnitude', 'phasediff'
    :return: string containing the BIDS modality folder
    """
    category_dir = {
        'T1w': 'anat',
        'T2w': 'anat',
        'FLAIR': 'anat',
        'bold': 'func',
        'dwi': 'dwi',
        'magnitude': 'fmap',
        'phasediff': 'fmap'}
    return category_dir[type]


def get_files_ext(type):
    """Give the BIDS file extensions to be fetched in function of the modality type

    :param type: the modality type amoung the following list: 'T1w', 'T2w', 'FLAIR', 'bold', 'dwi', 'magnitude', 'phasediff'
    :return: table of strings containing the extension of the expected files
    """
    files_ext = {
        'T1w': ['nii.gz', 'json'],
        'T2w': ['nii.gz', 'json'],
        'FLAIR': ['nii.gz', 'json'],
        'bold': ['nii.gz', 'json'],
        'dwi': ['nii.gz', 'json', 'bvals', 'bvecs'],
        'magnitude': ['nii.gz', 'json'],
        'phasediff': ['nii.gz', 'json']}
    return files_ext[type]
