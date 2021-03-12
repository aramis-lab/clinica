# coding: utf8

"""
This file contains a set of functional tests designed to check the correct execution of the pipeline and the
different functions available in Clinica
"""


import warnings
import sys
from testing_tools import create_list_hashes, compare_folders_with_hashes, compare_folders_structures
from testing_tools import clean_folder, compare_folders
from testing_tools import identical_subject_list, same_missing_modality_tsv
from os import pardir

# Determine location for working_directory
warnings.filterwarnings("ignore")


def test_run_Nifd2Bids(cmdopt):
    from clinica.iotools.converters.nifd_to_bids.nifd_to_bids import convert_clinical_data, convert_images
    from os.path import dirname, join, abspath
    import shutil

    root = join(dirname(abspath(__file__)), pardir, 'data', 'Nifd2Bids')

    clean_folder(join(root, 'out', 'bids'), recreate=True)
    clean_folder(join(root, 'out', 'clinical_data'), recreate=False)

    shutil.copytree(join(root, 'in', 'clinical_data'), join(root, 'out', 'clinical_data'))

    # Data location
    dataset_directory = join(root, 'in', 'unorganized')
    bids_directory = join(root, 'out', 'bids')
    clinical_data_directory = join(root, 'out', 'clinical_data')

    # Conversion
    to_convert = convert_images(dataset_directory, bids_directory, clinical_data_directory)
    convert_clinical_data(bids_directory, clinical_data_directory, to_convert)

    compare_folders_structures(bids_directory, join(root, 'ref', 'hashes_nifd.p'))

    clean_folder(join(root, 'out', 'bids'), recreate=True)
    clean_folder(join(root, 'out', 'clinical_data'), recreate=False)


def test_run_Oasis2Bids(cmdopt):
    from clinica.iotools.converters.oasis_to_bids.oasis_to_bids import OasisToBids
    from os.path import dirname, join, abspath

    root = join(dirname(abspath(__file__)), pardir, 'data', 'Oasis2Bids')

    clean_folder(join(root, 'out', 'bids'), recreate=True)

    # Data location
    dataset_directory = join(root, 'in', 'unorganized')
    bids_directory = join(root, 'out', 'bids')
    clinical_data_directory = join(root, 'in', 'clinical_data')

    oasis_to_bids = OasisToBids()
    oasis_to_bids.convert_images(dataset_directory, bids_directory)
    oasis_to_bids.convert_clinical_data(clinical_data_directory, bids_directory)

    compare_folders(join(root, 'out'), join(root, 'ref'),
                    shared_folder_name='bids')
    clean_folder(join(root, 'out', 'bids'), recreate=True)

def test_run_Oasis3ToBids(cmdopt):
    from clinica.iotools.converters.oasis3_to_bids.oasis3_to_bids import Oasis3ToBids
    from os.path import dirname, join, abspath

    root = join(dirname(abspath(__file__)), pardir, 'data', 'Oasis3ToBids')

    clean_folder(join(root, 'out', 'bids'), recreate=True)

    # Data location
    dataset_directory = join(root, 'in', 'unorganized')
    bids_directory = join(root, 'out', 'bids')
    clinical_data_directory = join(root, 'in', 'clinical_data')

    oasis_to_bids = Oasis3ToBids()
    oasis_to_bids.convert_images(dataset_directory, bids_directory)
    # oasis_to_bids.convert_clinical_data(clinical_data_directory, bids_directory)
    #
    # compare_folders(join(root, 'out'), join(root, 'ref'),
    #                 shared_folder_name='bids')
    # clean_folder(join(root, 'out', 'bids'), recreate=True)

def test_run_Adni2Bids(cmdopt):
    from clinica.iotools.converters.adni_to_bids.adni_to_bids import AdniToBids
    from os.path import dirname, join, abspath

    root = join(dirname(abspath(__file__)), pardir, 'data', 'Adni2Bids')

    clean_folder(join(root, 'out', 'bids'), recreate=True)

    adni_to_bids = AdniToBids()
    adni_to_bids.check_adni_dependencies()

    dataset_directory = join(root, 'in', 'unorganized_data')
    clinical_data_directory = join(root, 'in', 'clinical_data_25-04-19')
    bids_directory = join(root, 'out', 'bids')
    subjects_list = join(root, 'in', 'subjects.txt')
    modalities = ['T1', 'PET_FDG', 'PET_AMYLOID', 'PET_TAU', 'DWI', 'FLAIR', 'fMRI']
    adni_to_bids.convert_images(dataset_directory,
                                clinical_data_directory,
                                bids_directory,
                                subjects_list,
                                modalities)
    adni_to_bids.convert_clinical_data(clinical_data_directory, bids_directory)
    # Generate tree of output files
    compare_folders(join(root, 'out'), join(root, 'ref'),
                    shared_folder_name='bids')
    clean_folder(join(root, 'out', 'bids'), recreate=True)


def test_run_CreateSubjectSessionList(cmdopt):
    from os.path import join, dirname, abspath
    from os import remove
    from clinica.iotools.utils import data_handling as dt

    root = join(dirname(abspath(__file__)), pardir, 'data', 'CreateSubjectSessionList')

    # Set variables
    bids_directory = join(root, 'in', 'bids')
    output_directory = join(root, 'out')
    tsv_name = 'subject_session_list.tsv'

    # Create subject_session file
    dt.create_subs_sess_list(bids_directory, output_directory, tsv_name)

    # Comparison bitwise
    out_tsv = join(output_directory, tsv_name)
    ref_tsv = join(root, 'ref', tsv_name)
    assert identical_subject_list(out_tsv, ref_tsv)
    remove(out_tsv)


def test_run_CreateMergeFile(cmdopt):
    from os.path import join, dirname, abspath
    from os import remove
    from filecmp import cmp
    import shutil
    from clinica.iotools.utils import data_handling as dt

    root = join(dirname(abspath(__file__)), pardir, 'data', 'CreateMergeFile')

    bids_directory = join(root, 'in', 'bids')
    out_tsv = join(root, 'out', 'output_file.tsv')
    subject_session_tsv = join(root, 'in', 'subjects_sessions.tsv')

    clean_folder(join(root, 'out', 'caps'), recreate=False)
    shutil.copytree(join(root, 'in', 'caps'), join(root, 'out', 'caps'))
    caps_directory = join(root, 'out', 'caps')

    dt.create_merge_file(bids_directory,
                         out_tsv,
                         caps_dir=caps_directory,
                         pipelines=None,
                         atlas_selection=None,
                         pvc_restriction=None,
                         tsv_file=subject_session_tsv,
                         group_selection=None)
    # Comparison step
    ref_tsv = join(root, 'ref', 'output_file.tsv')
    assert cmp(out_tsv, ref_tsv)
    remove(out_tsv)


def test_run_ComputeMissingModalities(cmdopt):
    from os.path import join, dirname, abspath, exists
    from os import remove
    from clinica.iotools.utils import data_handling as dt

    root = join(dirname(abspath(__file__)), pardir, 'data', 'ComputeMissingMod')

    bids_directory = join(root, 'in', 'bids')
    output_directory = join(root, 'out')
    output_name = 'missing_modalities'

    dt.compute_missing_mods(bids_directory, output_directory, output_name)

    filenames = ['missing_modalities_ses-M00.tsv',
                 'missing_modalities_ses-M03.tsv',
                 'missing_modalities_ses-M06.tsv',
                 'missing_modalities_ses-M12.tsv',
                 'missing_modalities_ses-M24.tsv',
                 'missing_modalities_ses-M48.tsv']
    for i in range(len(filenames)):
        outname = join(root, 'out', filenames[i])
        refname = join(root, 'ref', filenames[i])
        if not exists(outname):
            raise FileNotFoundError('A file called ' + outname + ' should have been generated, but it does not exists')
        assert same_missing_modality_tsv(outname, refname)
        remove(join(root, 'out', filenames[i]))


def test_run_Aibl2Bids(cmdopt):
    from clinica.iotools.converters.aibl_to_bids.aibl_to_bids import convert_clinical_data, convert_images
    from os.path import dirname, join, abspath

    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, 'data', 'Aibl2Bids')

    dataset_directory = join(root, 'in', 'unorganized_data')
    clinical_data_directory = join(root, 'in', 'Data_extract_3.2.5')
    bids_directory = join(root, 'out', 'bids')

    clean_folder(join(root, 'out', 'bids'), recreate=True)

    # Perform conversion of dataset
    convert_images(dataset_directory, clinical_data_directory, bids_directory)
    convert_clinical_data(bids_directory, clinical_data_directory)

    # Evaluate difference between ref and out
    compare_folders(join(root, 'out'), join(root, 'ref'),
                    shared_folder_name='bids')
    clean_folder(join(root, 'out', 'bids'), recreate=True)


def test_run_CenterNifti(cmdopt):
    from clinica.iotools.utils.data_handling import center_all_nifti
    from os.path import dirname, join, abspath

    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, 'data', 'CenterNifti')

    clean_folder(join(root, 'out', 'bids_centered'), recreate=True)

    bids_dir = join(root, 'in', 'bids')
    output_dir = join(root, 'out', 'bids_centered')

    all_modalities = ['t1w', 'pet', 'dwi', 'magnitude', 'bold', 'flair', 't2', 'phasediff']

    center_all_nifti(bids_dir, output_dir, all_modalities, center_all_files=True)
    hashes_out = create_list_hashes(output_dir, extensions_to_keep=('.nii.gz', '.nii'))
    hashes_ref = create_list_hashes(join(root, 'ref', 'bids_centered'), extensions_to_keep=('.nii.gz', '.nii'))

    if hashes_out != hashes_ref:
        raise RuntimeError('Hashes of nii* files are different between out and ref')
    clean_folder(join(root, 'out', 'bids_centered'), recreate=False)
