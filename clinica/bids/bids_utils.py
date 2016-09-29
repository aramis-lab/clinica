from os import path
import os
import logging
from glob import glob
import fileinput
from shutil import copy
import nibabel


def remove_rescan(list_path):
    """
    Remove all the folders containing the keyword 'rescan' from
    a given list of folders.

    Args:
        list_path (str): the list of the files to analize.

    Returns:
        The list of non rescanned folders.

    """
    # Extract all the folders without the substring 'rescan'
    noResc_lst = [s for s in list_path if 'rescan' not in s]
    if len(noResc_lst) != len(list_path):
         for r_file in list(set(list_path)-set(noResc_lst)):
             logging.warning('Rescan found '+r_file+' Ignored.')

    return noResc_lst


def choose_correction(dir, to_consider, mod):
    """
    Decides what is the best file type to choose depending on a priority list.

    The selection criteria is based on the input list toConsider that contains
    the files to consider in order of priority.

    Args:
        dir: the directory containing the files.
        to_consider: the list of files to consider in order of priority.
        mod: the modality

    Returns:
        The best file (according to the list to_consider) for the given modality or:
        -1 if if the modality is missing.
        0 if none of the desided corrections is available..

    """
    # extract all the files availabe for a certain modality
    correction_list = remove_rescan(glob(path.join(dir, '*' + mod + '*')))
    if len(correction_list)==0:
        return -1
    if len(correction_list) == 1:
        return correction_list[0].split(os.sep)[-1]
    else:
        for i in range(0, len(to_consider)):
            if any(to_consider[i] in c for c in correction_list):
                return to_consider[i]
        return 0


def get_bids_suff():
    bids_suff = {
        'T1': '_T1w',
        'T2': '_T2w',
        'Flair': '_FLAIR',
        'MapPh': '_phasediff',
        'Map': '_magnitude',
        'fMRI': '_bold',
        'dwi': '_dwi'
    }
    return bids_suff


def convert_T1(t1_path, output_path, t1_bids_name, log=None):
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    copy(t1_path, path.join(output_path, t1_bids_name + (get_bids_suff())['T1'] + '.nii.gz'))


def convert_flair(folder_input, folder_output, name):
    flair_lst = remove_rescan(glob(path.join(folder_input,'*T2FLAIR*')))
    if len(flair_lst) == 1:
        if not os.path.exists(folder_output):
            os.mkdir(folder_output)
        flair_path = glob(path.join(flair_lst[0], '*.nii.gz*'))[0]
        copy(flair_path, path.join(folder_output, name + (get_bids_suff())['Flair'] + '.nii.gz'))
    elif len(flair_lst) == 0:
            logging.info('No FLAIR found for ' + folder_input)
            return -1
    elif len(flair_lst)>1:
            print('Impossible decide the FLAIR to convert. Computation aborted.')
            raise


def convert_fmri(folder_input, folder_output, name):
    """
    Look for fmri file(s) inside the input folder and convert it(them) into BIDS specification.
    Args:
        folder_input: the folder containing the fmri to convert
        folder_output: the output folder
        name:

    Returns:
        -1 in case that no fmri file is found within the folder
    """
    fmri_lst = remove_rescan(glob(path.join(folder_input, '*fMRI*')))
    if len(fmri_lst) > 0:
        os.mkdir(folder_output)
        fmri_file_path = glob(path.join(fmri_lst[0], '*.nii*'))[0]
        copy(fmri_file_path , path.join(folder_output, name + '_task-rest' + (get_bids_suff())['fMRI'] + '.nii.gz'))
    else:
        logging.info('Non fMRI found for ' + folder_input)
        return -1


def merge_DTI(folder_input, folder_output, name):
    """
    Merge all the DTI files of a given subject.

    For the merging only DTI folders containing all .nii.gz, .bval and .bvec are considered,
    otherwise the folder is ignored.

    Args:
        folder_input (str) : the folder containing all the DTI to merge.
        folder_output (str) : the folder where store the merged file.
        name (str): name to give to the merged file.

    Returns:
        -1 if the input folder doesn't contain any DTI folder.
    """
    img = []
    bval = []
    bvec = []
    #merger = Merge()
    dti_list = remove_rescan(glob(path.join(folder_input, '*DTI*')))
    incomp_folders = []
    nr_dti = len(dti_list)
    if nr_dti == 0:
        return -1
    else:
        if not os.path.exists(folder_output):
            os.mkdir(folder_output)
        for folder in dti_list:
            if len(glob(path.join(folder,'*.bval'))) != 0 and len(glob(path.join(folder,'*.bvec'))) != 0:
                img.append(glob(path.join(folder,'*.nii*'))[0])
                bval.append(glob(path.join(folder,'*.bval'))[0])
                bvec.append(glob(path.join(folder,'*.bvec'))[0])
            else:
                incomp_folders.append(folder)

        # if it has been found at least a DTI folder complete with bvec, bval and nii.gz
        if len(img) > 0:
            #merge all the .nii.gz file with fslmerge
            os.system('fslmerge -t '+path.join(folder_output,name+'.nii.gz')+' '+" ".join(img))
            # merger.inputs.in_files = img
            # merger.inputs.dimension = 't'
            # merger.inputs.output_type = 'NIFTI_GZ'
            # merger.cmdline
            #merge all the .bval files
            fin = fileinput.input(bval)
            fout = open(path.join(folder_output,name+'.bval'), 'w')
            for line in fin:
                fout.write(line)
            #merge all the .bvec files
            fin = fileinput.input(bvec)
            fout = open(path.join(folder_output, name + '.bvec'), 'w')
            for line in fin:
                fout.write(line)

        if len(incomp_folders) > 0:
            return incomp_folders
        else:
            return None