# coding: utf8


def likeliness_measure(file1, file2, threshold1, threshold2, display=False):
    """
    Function that compares 2 Nifti inputs, with 2 different thresholds.

    Args:
        (string) file1: path to first nifti input
        (string) file2: path to second nifti to compare
        (tuple) threshold1: defines the first criteria to meet: threshold1[0] defines the relative
                            difference between 2 voxels to be considered different (ex: 1e-4). threshold[1] defines
                            the maximum proportion of voxels that can different for the test to be negative.
        (tuple) threshold2: defines the second criteria to meet.
        (bool) display: If set to True, will display a useful graph to determine optimal threshold for the
                        comparison

    Returns:
        (bool) True if file1 and file2 can be considered similar enough (meeting criterion expressed in threshold1
               and threshold2). False otherwise.

    """
    import os

    import matplotlib.pyplot as plt
    import nibabel as nib
    import numpy as np

    print(" ** comparing " + os.path.basename(file1) + " **")
    data1 = nib.load(file1).get_fdata(dtype="float32")
    data1[np.isnan(data1)] = 0

    data2 = nib.load(file2).get_fdata(dtype="float32")
    data2[np.isnan(data2)] = 0

    # Get mask where data are 0 in data1 and data2
    mask = (data1 == 0) & (data2 == 0)
    data1[mask] = 1
    data2[mask] = 1
    metric = (2 * np.abs(data1 - data2)) / (np.abs(data1) + np.abs(data2))
    metric_flattened = np.ndarray.flatten(metric)

    # Display fig
    if display:
        thresholds = np.logspace(-8, 0, 20)
        percents = np.array(
            [np.sum((metric_flattened > T)) / metric_flattened.size for T in thresholds]
        )
        fig, ax = plt.subplots()
        ax.semilogx(thresholds, percents)
        ax.grid()
        plt.xlabel("Threshold of relative difference")
        plt.ylabel("Proportion of different voxels")
        plt.show()

    mask_different_voxels_cond1 = metric_flattened > threshold1[0]
    mask_different_voxels_cond2 = metric_flattened > threshold2[0]
    return (
        np.sum(mask_different_voxels_cond1) / metric_flattened.size < threshold1[1]
    ) & (np.sum(mask_different_voxels_cond2) / metric_flattened.size < threshold2[1])


def similarity_measure(file1, file2, threshold):
    """
    Function that compares 2 Nifti inputs using a correlation metric. Nifti are equals if correlation gives

    Args:
        (string) file1: path to first nifti input
        (string) file2: path to second nifti to compare
        (float) threshold

    Returns:
        (bool) True if file1 and file2 can be considered similar enough. (superior than threshold)

    """
    import nipype.pipeline.engine as npe
    import numpy as np
    from nipype.algorithms.metrics import Similarity

    # Node similarity (nipy required)
    img_similarity = npe.Node(name="img_similarity", interface=Similarity())
    img_similarity.inputs.metric = "cc"  # stands for correlation coefficient
    img_similarity.inputs.volume1 = file1
    img_similarity.inputs.volume2 = file2
    res = img_similarity.run()

    return np.mean(res.outputs.similarity) > threshold


def identical_subject_list(sub_ses_list1, sub_ses_list2):
    """
    Function that ensures that both subject_session files are describing the same list

    Args:
        (string) sub_ses_list1: path to first nifti input
        (string) sub_ses_list2: path to second nifti to compare

    Returns:
        (bool) True if sub_ses_list1 and sub_ses_list2 contains the same sessions

    """

    def is_included(list1, list2):
        from pandas import read_csv

        # Read csv files
        readlist1 = read_csv(list1, sep="\t")
        readlist2 = read_csv(list2, sep="\t")

        # If columns are different, files are different
        if list(readlist1.columns) != list(readlist2.columns):
            return False
        else:

            # Extract subject and corresponding session names
            subjects1 = list(readlist1.participant_id)
            sessions1 = list(readlist1.session_id)
            subjects2 = list(readlist2.participant_id)
            sessions2 = list(readlist2.session_id)

            if len(subjects1) != len(subjects2):
                return False
            else:
                for i in range(len(subjects1)):
                    current_sub = subjects1[i]
                    current_ses = sessions1[i]
                    # Compute all the indices in the second list corresponding to the current subject
                    idx_same_sub = [
                        j for j in range(len(subjects2)) if subjects2[j] == current_sub
                    ]
                    if len(idx_same_sub) == 0:  # Current subject not found in
                        return False
                    ses_same_sub = [sessions2[idx] for idx in idx_same_sub]
                    if current_ses not in ses_same_sub:
                        return False
        return True

    # The operation is performed both sides because is_included(list1, list2) != is_included(list2, list1)
    return is_included(sub_ses_list1, sub_ses_list2) & is_included(
        sub_ses_list2, sub_ses_list1
    )


def same_missing_modality_tsv(file1, file2):
    """
    Function that is used to compare 2 TSV files generated by the iotool ComputeMissingModalities.

    Only fields participant_id, pet, t1w, func_task - rest are compared. Line order does not matter.


    Args:
        (string) file1: path to first tsv
        (string) file2: path to second tsv

    Returns:
        (bool) True if file1 and file2 contains the same information
    """
    import pandas as pds

    # Read dataframe with pandas
    df1 = pds.read_csv(file1, sep="\t")
    df2 = pds.read_csv(file2, sep="\t")

    # Extract data and form lists for both files
    subjects1 = list(df1.participant_id)
    pet_AV45_1 = list(df1["pet_acq-AV45"])
    pet_FDG_1 = list(df1["pet_acq-FDG"])
    t1w1 = list(df1.t1w)
    func_task_rest1 = list(df1["func_task-rest"])

    subjects2 = list(df2.participant_id)
    pet_AV45_2 = list(df2["pet_acq-AV45"])
    pet_FDG_2 = list(df2["pet_acq-FDG"])
    t1w2 = list(df2.t1w)
    func_task_rest2 = list(df2["func_task-rest"])

    # Subjects are sorted in alphabetical order. The same permutation of element is applied on each column
    subjects1_sorted, pet_AV45_1 = (
        list(t) for t in zip(*sorted(zip(subjects1, pet_AV45_1)))
    )
    subjects2_sorted, pet_AV45_2 = (
        list(t) for t in zip(*sorted(zip(subjects2, pet_AV45_2)))
    )
    subjects1_sorted, pet_FDG_1 = (
        list(t) for t in zip(*sorted(zip(subjects1, pet_FDG_1)))
    )
    subjects2_sorted, pet_FDG_2 = (
        list(t) for t in zip(*sorted(zip(subjects2, pet_FDG_2)))
    )
    subjects1_sorted, t1w1 = (list(t) for t in zip(*sorted(zip(subjects1, t1w1))))
    subjects2_sorted, t1w2 = (list(t) for t in zip(*sorted(zip(subjects2, t1w2))))
    subjects1_sorted, func_task_rest1 = (
        list(t) for t in zip(*sorted(zip(subjects1, func_task_rest1)))
    )
    subjects2_sorted, func_task_rest2 = (
        list(t) for t in zip(*sorted(zip(subjects2, func_task_rest2)))
    )

    # Test is positive when all the sorted list s are equals
    return (
        (subjects1_sorted == subjects2_sorted)
        & (pet_AV45_1 == pet_AV45_2)
        & (pet_FDG_1 == pet_FDG_2)
        & (t1w1 == t1w2)
        & (func_task_rest1 == func_task_rest2)
    )


def compare_folders(out, ref, shared_folder_name):
    from filecmp import cmp
    from os import remove
    from os.path import join

    out_txt = join(out, "out_folder.txt")
    ref_txt = join(ref, "ref_folder.txt")

    list_files(join(out, shared_folder_name), filename=out_txt)
    list_files(join(ref, shared_folder_name), filename=ref_txt)

    # Compare them
    if not cmp(out_txt, ref_txt):
        with open(out_txt, "r") as fin:
            out_message = fin.read()
        with open(ref_txt, "r") as fin:
            ref_message = fin.read()
        remove(out_txt)
        remove(ref_txt)
        raise ValueError(
            "Comparison of out and ref directories shows mismatch :\n "
            "OUT :\n" + out_message + "\n REF :\n" + ref_message
        )

    # Clean folders
    remove(out_txt)
    remove(ref_txt)


def list_files(startpath, filename=None):
    """

    Args:
        startpath: starting point for the tree listing. Does not list hidden
        files (to avoid problems with .DS_store for example
        filename: if None, display to stdout, otherwise write in the file

    Returns:
        void
    """
    from os import remove, sep, walk
    from os.path import abspath, basename, exists, expanduser, expandvars

    if exists(filename):
        remove(filename)

    expanded_path = abspath(expanduser(expandvars(startpath)))
    for root, dirs, files in walk(expanded_path):
        level = root.replace(startpath, "").count(sep)
        indent = " " * 4 * (level)
        rootstring = "{}{}/".format(indent, basename(root))
        # Do not deal with hidden files
        if not basename(root).startswith("."):
            if filename is not None:
                # 'a' stands for 'append' rather than 'w' for 'write'. We must
                # manually jump line with \n otherwise everything is
                # concatenated
                with open(filename, "a") as fin:
                    fin.write(rootstring + "\n")
            else:
                print(rootstring)
            subindent = " " * 4 * (level + 1)
            for f in files:
                filestring = "{}{}".format(subindent, f)
                if not basename(f).startswith("."):
                    if filename is not None:
                        with open(filename, "a") as fin:
                            fin.write(filestring + "\n")
                    else:
                        print(filestring)


def clean_folder(path, recreate=True):
    from os import makedirs
    from os.path import abspath, exists
    from shutil import rmtree

    abs_path = abspath(path)
    if exists(abs_path):
        rmtree(abs_path)
    if recreate:
        makedirs(abs_path)


def clean_PETLinear(caps_path):
    from os import listdir, makedirs, path
    from os.path import abspath, exists
    from shutil import rmtree

    # Handle subjects subdirectory
    abs_path_subjs = path.join(caps_path, "subjects")

    # Iterate over the subject directories
    for subj in listdir(abs_path_subjs):
        subj_dir = path.join(abs_path_subjs, subj)
        if path.isdir(subj_dir):
            # Iterate over session directories
            for sess in listdir(subj_dir):
                sess_dir = path.join(subj_dir, sess)
                if path.isdir(sess_dir):

                    path_pet_linear = path.join(sess_dir, "pet_linear")
                    clean_folder(path_pet_linear, recreate=False)

    # TODO: handle groups subdirectory


def create_list_hashes(path_folder, extensions_to_keep=(".nii.gz", ".tsv", ".json")):
    """
    Computes a dictionary of files with their corresponding hashes

        Args:
            (string) path_folder: starting point for the tree listing.
            (tuple) extensions_to_keep: files with these extensions will have their hashes computed and tracked

        Returns:
            (dictionary) all_files: a dictionary of the form {/path/to/file.extension: hash(file.extension)}
    """
    import hashlib
    import os

    def file_as_bytes(file):
        with file:
            return file.read()

    all_files = []
    for subdir, dirs, files in os.walk(path_folder):
        files.sort()
        for file in files:
            if file.lower().endswith(extensions_to_keep):
                all_files.append(os.path.join(subdir, file))

    dict_hashes = {
        fname[len(path_folder) :]: str(
            hashlib.md5(file_as_bytes(open(fname, "rb"))).digest()
        )
        for fname in all_files
    }
    return dict_hashes


def compare_folders_with_hashes(path_folder, list_hashes):
    """
    Compares the files of a folder against a reference

        Args:
            (string) path_folder: starting point for the tree listing.
            (dictionary) list_hashes: a dictionary of the form {/path/to/file.extension: hash(file.extension)}
    """
    import pickle

    hashes_check = pickle.load(open(list_hashes, "rb"))
    hashes_new = create_list_hashes(path_folder)

    if hashes_check != hashes_new:
        error_message1 = ""
        error_message2 = ""
        for key in hashes_check:
            if key not in hashes_new:
                error_message1 += "{0} not found !\n".format(key)
            elif hashes_check[key] != hashes_new[key]:
                error_message2 += "{0} does not match the reference file !\n".format(
                    key
                )
        raise ValueError(error_message1 + error_message2)


def compare_folders_structures(path_folder, list_hashes):
    """
    Compares the structure of a folder against a reference

        Args:
            (string) path_folder: starting point for the tree listing.
            (dictionary) list_hashes: a dictionary of the form {/path/to/file.extension: hash(file.extension)}
    """
    import pickle

    hashes_check = pickle.load(open(list_hashes, "rb"))
    hashes_new = create_list_hashes(path_folder)

    if list(hashes_check).sort() != list(hashes_new).sort():
        error_message1 = ""
        error_message2 = ""
        for key in hashes_check:
            if key not in hashes_new:
                error_message1 += "{0} not found !\n".format(key)
        for key in hashes_new:
            if key not in hashes_check:
                error_message2 += "{0}'s creation was not expected !\n".format(key)

        raise ValueError(error_message1 + error_message2)
