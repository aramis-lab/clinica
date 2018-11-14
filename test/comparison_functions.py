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
    import nibabel as nib
    import numpy as np
    import os
    import matplotlib.pyplot as plt

    print(' ** comparing ' + os.path.basename(file1) + ' **')
    data1 = nib.load(file1).get_data()
    data1[data1 != data1] = 0

    data2 = nib.load(file2).get_data()
    data2[data2 != data2] = 0

    # Get mask where data are 0 in data1 and data2
    mask = (data1 == 0) & (data2 == 0)
    data1[mask] = 1
    data2[mask] = 1
    metric = (2 * np.abs(data1 - data2)) / (np.abs(data1) + np.abs(data1))
    metric_flattened = np.ndarray.flatten(metric)
    thresholds = np.logspace(-8, 0, 20)
    percents = np.array([np.sum((metric_flattened > T)) / metric_flattened.size for T in thresholds])

    # Display fig
    if display:
        fig, ax = plt.subplots()
        ax.semilogx(thresholds, percents)
        ax.grid()
        plt.xlabel('Threshold of relative difference')
        plt.ylabel('Proportion of different voxels')
        plt.show()

    mask_different_voxels_cond1 = metric_flattened > threshold1[0]
    mask_different_voxels_cond2 = metric_flattened > threshold2[0]
    return (np.sum(mask_different_voxels_cond1) / metric_flattened.size < threshold1[1]) & \
           (np.sum(mask_different_voxels_cond2) / metric_flattened.size < threshold2[1])


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
    import numpy as np
    import nipype.pipeline.engine as npe
    from nipype.algorithms.metrics import Similarity

    # Node similarity (nipy required)
    img_similarity = npe.Node(name='img_similarity',
                              interface=Similarity())
    img_similarity.inputs.metric = 'cc'  # stands for correlation coefficient
    img_similarity.inputs.volume1 = file1
    img_similarity.inputs.volume2 = file2
    res = img_similarity.run()

    return np.mean(res.outputs.similarity) > threshold