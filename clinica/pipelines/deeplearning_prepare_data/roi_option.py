

def extract_roi(input_path, basedir, masks_location, input_pattern, roi_list, uncrop_output):
    """Extracts regions of interest defined by masks

    This function extracts regions of interest from the preprocessed nifti image.
    The regions are defined using binary masks that must be located in the CAPS
    at `masks/roi_based/t1_linear`.

    Args:
        input_path: path to the tensor version of the nifti MRI.
        basedir: path to the extracted files.
        roi_list: list of the names of the regions that will be extracted.
        uncrop_output: if True, the final region is not cropped.

    Returns:
        file: multiple tensors saved on the disk, suffixes corresponds to
            indexes of the patches. Same location than input file.
    """
    import torch
    import nibabel as nib
    import os
    import numpy as np
    from os import path

    crop_desc = ""
    if not uncrop_output:
        crop_desc = "_desc-Crop"

    image_tensor = torch.load(input_path)

    input_tensor_filename = os.path.basename(input_path)
    txt_idx = input_tensor_filename.rfind("_")
    it_filename_prefix = input_tensor_filename[0:txt_idx]
    it_filename_suffix = input_tensor_filename[txt_idx:]

    output_roi = []
    for index_roi, roi in enumerate(roi_list):
        roi_path = path.join(masks_location, "roi-%s_%s_mask.nii.gz" % (input_pattern, roi))
        roi_mask = nib.load(roi_path).get_data()
        if len(roi_mask.shape) == 4:
            roi_mask = roi_mask[0]

        extracted_roi = image_tensor * roi_mask
        if not uncrop_output:
            extracted_roi = extracted_roi[np.ix_(roi_mask.any((1, 2)), roi_mask.any((0, 2)), roi_mask.any((0, 1)))]
        # save into .pt format
        output_roi.append(
            os.path.join(
                basedir,
                "%sroi-%s%s%s" % (it_filename_prefix, roi, crop_desc, it_filename_suffix)
            )
        )
        torch.save(extracted_roi.clone(), output_roi[index_roi])

    return output_roi


def check_mask_list(masks_location, roi_list):
    from os import path

    for roi in enumerate(roi_list):
        roi_path = path.join(masks_location, "roi-%s_mask.nii.gz" % roi)
        if not path.exists(roi_path):
            raise ValueError('The ROI wanted do not correspond to a mask in the CAPS directory.'
                             'It should be written at the following path: %s' % roi_path)


if __name__ == "__main__":
    import argparse
    from os import path
    import os

    from ...utils.inputs import clinica_file_reader

    parser = argparse.ArgumentParser(description="Temporary parser for the extraction of ROI tensors.")

    parser.add_argument("caps_directory", type=str, help="path to the CAPS directory.")
    parser.add_argument("--roi_list", type=str, nargs="+", default=None,
                        help="List of regions to be extracted")
    parser.add_argument("--uncrop_output", action="store_true", default=False,
                        help="Disable cropping option so the output tensors have the same size than the whole image.")
    args = parser.parse_args()

    preprocessing = "t1_linear"
    input_pattern = "space-MNI152NLin2009cSym_desc-Crop_res-1x1x1"

    # Load the corresponding masks
    masks_location = path.join(args.caps_directory, 'masks', 'roi_based', preprocessing)

    if args.roi_list is None:
        args.roi_list = [mask.split("_")[0].split("-")[1] for mask in os.listdir(masks_location)]
        if len(args.roi_list) == 0:
            raise ValueError('No ROI was specified, and no mask was found at %s' % masks_location)
    else:
        check_mask_list(masks_location, args.roi_list)

    for subject in os.listdir(path.join(args.caps_directory, "subjects")):
        subject_path = path.join(args.caps_directory, "subjects", subject)
        for session in os.listdir(subject_path):
            input_path = clinica_file_reader([subject], [session], args.caps_directory,
                                             {"needed_pipeline": preprocessing,
                                              "pattern": input_pattern})[0]
            output_path = path.join(subject_path, session, "deeplearning_prepare_data", preprocessing)
            extract_roi(input_path, output_path, masks_location, input_pattern,
                        roi_list=args.roi_list, uncrop_output=args.uncrop_output)



