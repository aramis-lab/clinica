# coding: utf8

"""
This file contains a set of functional tests designed to check the correct execution of the pipeline and the
different functions available in Clinica
"""

__author__ = "Arnaud Marcoux"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Arnaud Marcoux"]
__license__ = "See LICENSE.txt file"
__version__ = "0.2.0"
__maintainer__ = "Arnaud Marcoux"
__email__ = "arnaud.marcoux@inria.fr"
__status__ = "Development"


import warnings
import sys
from .testing_tools import clean_folder, compare_folders
from .testing_tools import identical_subject_list, same_missing_modality_tsv
from os import pardir

# Determine location for working_directory
warnings.filterwarnings("ignore")

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


def test_run_Adni2Bids(cmdopt):
    from clinica.iotools.converters.adni_to_bids.adni_to_bids import AdniToBids
    from os.path import dirname, join, abspath

    root = join(dirname(abspath(__file__)), pardir, 'data', 'Adni2Bids')

    clean_folder(join(root, 'out', 'bids'), recreate=True)

    adni_to_bids = AdniToBids()
    adni_to_bids.check_adni_dependencies()

    dataset_directory = join(root, 'in', 'unorganized_data')
    clinical_data_directory = join(root, 'in', 'ADNI_clinical_data_new')
    bids_directory = join(root, 'out', 'bids')
    subjects_list = join(root, 'in', 'subjects.txt')
    # modalities = ['T1', 'PET_FDG', 'PET_AV45', 'DWI', 'fMRI']
    modalities = ['T1', 'PET_FDG', 'PET_AV45', 'DWI', 'FLAIR']
    adni_to_bids.convert_images(dataset_directory,
                                clinical_data_directory,
                                bids_directory,
                                subjects_list,
                                modalities)

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

