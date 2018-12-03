# coding: utf8


def permutation_test(vector_1, vector_2, number_of_permutation, tails=2):
    """
    This function take two vectors as entry and return the p_value
    of a permutation test.
    If permutation not done, p_value is set to 1.
    This function performs a two tails permutation on *non-nul values*
    If tails is set to 1 the contrast is supposed to be vector_1 > vector_2

    WARNING: here 0 are not removed from the dataset
    If tails is set to 1 the contrast is supposed to be vector_1 > vector_2
    """
    import numpy as np

    def exact_mc_perm_test(xs, ys, nmc, tails=2):
        n, k = len(xs), 0.
        if tails == 2:
            diff = np.abs(np.mean(xs) - np.mean(ys))
        if tails == 1:
            diff = np.mean(xs) - np.mean(ys)
        zs = np.concatenate([xs, ys])
        for j in range(nmc):
            np.random.shuffle(zs)
            if tails == 2:
                diff_random = np.abs(np.mean(zs[:n]) - np.mean(zs[n:]))
            if tails == 1:
                diff_random = np.mean(zs[:n]) - np.mean(zs[n:])
            k += diff < diff_random
        return k/float(nmc)

    p_value = 1.

    x = vector_1
    y = vector_2

    if x != np.array([]) and y != np.array([]):
        if x.sum() != 0 or y.sum() != 0:
            p_value = exact_mc_perm_test(x, y, number_of_permutation, tails)

    return p_value


def t_test(vector_1, vector_2, tails=2):
    """
    This function takes two vectors as entry and return the p-value
    of a t_test.
    If the t_test not done (only zeros in both vectors), p_value is set to 1.
    This function performs a t_test on *non-nul values*.
    If tails is set to 1 the contrast is supposed to be vector_1 > vector_2.

    WARNING: here 0 are not removed from the dataset
    If tails is set to 1 the contrast is supposed to be vector_1 > vector_2
    """

    import numpy as np
    from scipy.stats import ttest_ind

    p_value = 1.

    X = vector_1
    Y = vector_2

    if X != np.array([]) and Y != np.array([]):
        p_value = ttest_ind(X, Y, equal_var=False).pvalue
        if tails == 1:
            x = np.array(X)
            array_y = np.array(Y)
            if x.mean() > array_y.mean():
                p_value = float(p_value)/2.0
            else:
                p_value = 1 - float(p_value)/2.0

    return p_value


def mann_whitney(vector_1, vector_2, tails=2):
    """
    This function take two vectors as entry and return
    the p_value of a  Mann_Whitney U test.
    If the U test not done (only 0 in both vectors), p_value is set to 1.
    This function performs a U test on *non-nul values*

    WARNING: here 0 are not removed from the dataset
    If tails is set to 1 the contrast is supposed to be vector_1 > vector_2
    """

    import numpy as np
    from scipy.stats import mannwhitneyu

    p_value = 1.

    x = vector_1
    non_zero_x = [x[k] for k in np.nonzero(x)[0].tolist()]

    y = vector_2
    non_zero_y = [y[k] for k in np.nonzero(y)[0].tolist()]

    if x != np.array([]) and y != np.array([]):
        if non_zero_x != np.array([]) or non_zero_y != np.array([]):
            p_value = mannwhitneyu(x, y).pvalue
            if tails == 2:
                p_value = 2*p_value

    return p_value


def fdr_correction_matrix(p_value_matrix, template=None):
    """
    This function take a p value matrix as entry and return the corrected
    p_value for False Rate Discovery.
    If not all statistical tests have been performed (typically in DTI at
    a absent connection) a template matrix (which is a binary metrix with 1
    if the test is performed and 0 else) with the same shape as p_value_matrix
    input of the actually performed test can be provide at input.
    """

    import numpy as np
    from mne.stats import fdr_correction

    if type(template) == type(p_value_matrix):
        if p_value_matrix.shape != template.shape:
            raise IOError(
                'p_value_matrix and template should have the same shape.')

    if type(template) == type(p_value_matrix):
        p_value_corrected = np.ones(p_value_matrix.shape)
        reject_test = np.zeros(p_value_matrix.shape, dtype=bool)
        eff_p_value = []
        index_of_eff_p_value = []
        for i in np.arange(0, p_value_matrix.shape[0]):
            for j in np.arange(0, i):
                if template[j, i] == 1:
                    eff_p_value += [p_value_matrix[j, i]]
                    index_of_eff_p_value += [(j, i)]
        reject, p_corrected = fdr_correction(eff_p_value)
        for i, corrected in enumerate(p_corrected):
            p_value_corrected[
                index_of_eff_p_value[i][0], index_of_eff_p_value[i][1]
            ] = corrected
            reject_test[
                index_of_eff_p_value[i][0], index_of_eff_p_value[i][1]
            ] = reject[i]
    elif not template:
        reject_test, p_value_corrected = fdr_correction(p_value_matrix)
    else:
        raise IOError('template input should be an numpy array or None.')
    return reject_test, p_value_corrected


def create_new_feature_tsv(subjects_visits_tsv, bids_dir, dest_tsv, added_features):
    """
        This func is to add new features(columns) from the subjects_visits_list
        TSV file, and use the generated file in the statistical analysis
        added_features : list of str, the added features that you want to write
        into the tsv file, e.g. ['age_bl', 'sex']

    Args:
        subjects_visits_tsv: tsv files containing just participant_id and session_id columns
        bids_dir: BIDS directory
        dest_tsv: the destination tsv file
        added_features: a list containing all the new columns from the participant.tsv in BIDS directory.

    Returns:

    """
    import pandas as pd
    from os.path import join, isfile
    import logging
    from pandas import concat
    logging.basicConfig(level=logging.DEBUG)

    if not isfile(join(bids_dir, 'participants.tsv')):
        raise Exception('participants.tsv not found')
    sub_set = pd.io.parsers.read_csv(subjects_visits_tsv, sep='\t')
    all_set = pd.io.parsers.read_csv(join(bids_dir, 'participants.tsv'), sep='\t')
    all_set.set_index("participant_id", inplace=True)
    selected_subj = all_set.loc[list(sub_set.participant_id)]
    if sub_set.shape[0] != selected_subj.shape[0]:
        missing_subj = set(list(sub_set.participant_id)) - set(list(selected_subj.participant_id))
        msg = "The missing subjects are %s" % list(missing_subj)
        logging.info(msg)
        raise Exception('There are subjects which are not included in full dataset, please check it out')

    new_features = selected_subj[added_features]
    new_features.reset_index(inplace=True, drop=True)
    all_features = concat([sub_set, new_features], axis=1)
    all_features.to_csv(dest_tsv, sep='\t', index=False, encoding='utf-8')


def statistics_on_atlas(in_normalized_map, in_atlas, out_file=None):
    """
    Compute statistics of a map on an atlas.

    Given an atlas image with a set of ROIs, this function computes the mean of
    a normalized map (e.g. GM segmentation, FA map from DTI, etc.) on each ROI.

    Args:
        in_normalized_map (str): File containing a scalar image registered
            on the atlas.
        in_atlas (:obj: AbstractClass): An atlas with a set of ROI. These ROI
            are used to compute statistics.
        out_file (Optional[str]): Name of the output file.

    Returns:
        out_file (str): TSV file containing the statistics (content of the
            columns: label, mean scalar, std of the scalar', number of voxels).
    """
    from clinica.utils.atlas import AtlasAbstract
    import nibabel as nib
    import numpy as np
    import pandas
    import os.path as op
    from clinica.utils.stream import cprint

    if not isinstance(in_atlas, AtlasAbstract):
        raise Exception("Atlas element must be an AtlasAbstract type")

    if out_file is None:
        fname, ext = op.splitext(op.basename(in_normalized_map))
        if ext == ".gz":
            fname, ext2 = op.splitext(fname)
            ext = ext2 + ext
        out_file = op.abspath("%s_statistics_%s.tsv"
                              % (fname, in_atlas.get_name_atlas()))

    atlas_labels = nib.load(in_atlas.get_atlas_labels())
    atlas_labels_data = atlas_labels.get_data()

    img = nib.load(in_normalized_map)
    img_data = img.get_data()

    atlas_correspondance = pandas.io.parsers.read_csv(in_atlas.get_tsv_roi(), sep='\t')
    label_name = list(atlas_correspondance.roi_name)
    label_value = list(atlas_correspondance.roi_value)  # TODO create roi_value column in lut_*.txt and remove irrelevant RGB information

    mean_signal_value = []
    for label in label_value:
        current_mask_label = atlas_labels_data == label
        masked_data = np.array(img_data, copy=True)
        masked_data[np.invert(current_mask_label)] = 0
        mean_signal_value.append(np.sum(masked_data) / np.sum(current_mask_label))

    try:
        data = pandas.DataFrame({'label_name': label_name,
                                 'mean_scalar': mean_signal_value
                                 })
        data.to_csv(out_file, sep='\t', index=True, encoding='utf-8')
    except Exception as e:
        cprint("Impossible to save %s with pandas" % out_file)
        raise e

    return out_file
