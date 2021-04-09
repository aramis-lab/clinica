templates_dict = {"t1-linear": "MNI152NLin2009cSym", "t1-extensive": "Ixi549Space"}


def extract_roi(
    input_path,
    basedir,
    masks_location,
    mask_pattern,
    cropped_input,
    roi_list,
    uncrop_output,
):
    """Extracts regions of interest defined by masks

    This function extracts regions of interest from preprocessed nifti images.
    The regions are defined using binary masks that must be located in the CAPS
    at `masks/tpl-<template>`.

    Args:
        input_path: path to the tensor version of the nifti MRI.
        basedir: path to the extracted files.
        masks_location: path to the masks
        mask_pattern: pattern to identify the masks
        cropped_input: if the input is cropped or not (contains desc-Crop)
        roi_list: list of the names of the regions that will be extracted.
        uncrop_output: if True, the final region is not cropped.

    Returns:
        file: multiple tensors saved on the disk, suffixes corresponds to
            indexes of the patches. Same location than input file.
    """
    import os

    import nibabel as nib
    import numpy as np
    import torch

    image_array = nib.load(input_path).get_fdata()
    image_tensor = torch.from_numpy(image_array).unsqueeze(0).float()

    input_tensor_filename = os.path.basename(input_path)

    sub_ses_prefix = "_".join(input_tensor_filename.split("_")[0:3:])
    if not sub_ses_prefix.endswith("_T1w"):
        sub_ses_prefix = "_".join(input_tensor_filename.split("_")[0:2:])
    input_suffix = input_tensor_filename.split("_")[-1].split(".")[0]

    output_roi = []
    for index_roi, roi in enumerate(roi_list):
        mask_path = find_mask_path(masks_location, roi, mask_pattern, cropped_input)
        mask_np = nib.load(mask_path).get_fdata()
        if len(mask_np.shape) == 3:
            mask_np = mask_np[np.newaxis, :]

        extracted_roi = image_tensor * mask_np
        if not uncrop_output:
            extracted_roi = extracted_roi[
                np.ix_(
                    mask_np.any((1, 2, 3)),
                    mask_np.any((0, 2, 3)),
                    mask_np.any((0, 1, 3)),
                    mask_np.any((0, 1, 2)),
                )
            ]
        extracted_roi = extracted_roi.float()
        # save into .pt format
        output_pattern = compute_output_pattern(mask_path, not uncrop_output)
        output_roi.append(
            os.path.join(
                basedir, f"{sub_ses_prefix}_{output_pattern}_{input_suffix}.pt"
            )
        )
        os.makedirs(output_path, exist_ok=True)
        torch.save(extracted_roi.clone(), output_roi[index_roi])

    return output_roi


def check_mask_list(masks_location, roi_list, mask_pattern, cropping):
    import nibabel as nib
    import numpy as np

    for roi in roi_list:
        roi_path = find_mask_path(masks_location, roi, mask_pattern, cropping)
        if roi_path is None:
            raise ValueError(
                f"The ROI wanted do not correspond to a mask in the CAPS directory."
                f"The mask should include the following pattern: {mask_pattern}."
            )
        roi_mask = nib.load(roi_path).get_fdata()
        mask_values = set(np.unique(roi_mask))
        if mask_values != {0, 1}:
            raise ValueError(
                "The ROI masks used should be binary (composed of 0 and 1 only)."
            )


def find_mask_path(masks_location, roi, mask_pattern, cropping):
    """Finds masks corresponding to the pattern asked and containing the adequate cropping description"""
    from glob import glob
    from os import path

    candidates = glob(
        path.join(masks_location, f"*{mask_pattern}*roi-{roi}*_mask.nii.gz")
    )
    if cropping:
        candidates = [mask for mask in candidates if "_desc-Crop_" in mask]
    else:
        candidates = [mask for mask in candidates if "_desc-Crop_" not in mask]

    if len(candidates) == 0:
        return None
    else:
        return min(candidates, key=len)


def compute_output_pattern(mask_path, crop_output):
    """
    Computes the output pattern of the region cropped (without the source file prefix)

    Args:
        mask_path: path to the masks
        crop_output: If True the output is cropped, and the descriptor CropRoi must exist
    Returns:
        the output pattern
    """
    from os import path

    mask_filename = path.basename(mask_path)
    template_id = mask_filename.split("_")[0].split("-")[1]
    mask_descriptors = mask_filename.split("_")[1:-2:]
    roi_id = mask_filename.split("_")[-2].split("-")[1]
    if "desc-Crop" not in mask_descriptors and crop_output:
        mask_descriptors = ["desc-CropRoi"] + mask_descriptors
    elif "desc-Crop" in mask_descriptors:
        mask_descriptors = [
            descriptor for descriptor in mask_descriptors if descriptor != "desc-Crop"
        ]
        if crop_output:
            mask_descriptors = ["desc-CropRoi"] + mask_descriptors
        else:
            mask_descriptors = ["desc-CropImage"] + mask_descriptors

    mask_pattern = "_".join(mask_descriptors)
    output_pattern = f"space-{template_id}_{mask_pattern}_roi-{roi_id}"

    return output_pattern


if __name__ == "__main__":
    import argparse
    import os
    from os import path

    from clinica.utils.input_files import T1W_EXTENSIVE, T1W_LINEAR, T1W_LINEAR_CROPPED

    from ...utils.exceptions import ClinicaCAPSError
    from ...utils.inputs import clinica_file_reader

    parser = argparse.ArgumentParser(
        description="Temporary parser for the extraction of ROI tensors."
    )

    parser.add_argument("caps_directory", type=str, help="path to the CAPS directory.")
    parser.add_argument(
        "modality",
        help="""For which modality the tensor will be extracted.
                        't1-linear': images preprocessed with t1-linear pipeline.
                        't1-extensive': images preprocessed with t1-extensive pipeline.
                        'custom': find images with a custom suffix in their filename and
                        transform them to tensor format.""",
        choices=["t1-linear", "t1-extensive", "custom"],
        default="t1-linear",
    )
    parser.add_argument(
        "-uui",
        "--use_uncropped_image",
        help="""Use the uncropped image instead of the cropped image generated by t1-linear.""",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--roi_list",
        type=str,
        nargs="+",
        default=None,
        help="List of regions to be extracted",
    )
    parser.add_argument(
        "--uncrop_output",
        action="store_true",
        default=False,
        help="Disable cropping option so the output tensors have the same size than the whole image.",
    )
    parser.add_argument(
        "-cn",
        "--custom_suffix",
        help="""Custom suffix filename, e.g.:
            'graymatter_space-Ixi549Space_modulated-off_probability.nii.gz', or
            'segm-whitematter_probability.nii.gz'
            """,
        type=str,
        default="",
    )
    parser.add_argument(
        "-ct",
        "--custom_template",
        help="""Name of the template used when modality is set to custom.""",
        type=str,
        default=None,
    )
    parser.add_argument(
        "-cmp",
        "--custom_mask_pattern",
        help="""If given will select only the masks containing the string given.
                The mask with the shortest name is taken.
                This argument is taken into account only of the modality is custom.""",
        type=str,
        default=None,
    )
    args = parser.parse_args()

    if args.modality == "custom":
        template = args.custom_template
        if template is None:
            raise ValueError(
                "A custom template must be defined when the modality is set to custom."
            )
    else:
        template = templates_dict[args.modality]

    if args.modality == "t1-linear" and not args.use_uncropped_image:
        cropping = True
    else:
        cropping = False

    if args.modality == "t1-linear":
        if args.use_uncropped_image:
            FILE_TYPE = T1W_LINEAR
        else:
            FILE_TYPE = T1W_LINEAR_CROPPED
        output_folder = "t1_linear"
        args.mask_pattern = "res-1x1x1"
    elif args.modality == "t1-extensive":
        FILE_TYPE = T1W_EXTENSIVE
        output_folder = "t1_extensive"
        args.mask_pattern = "desc-skullstripped"
    else:
        FILE_TYPE = {
            "pattern": f"*{args.custom_suffix}",
            "description": "Custom suffix",
        }
        output_folder = "custom"

    # Load the corresponding masks
    masks_location = path.join(args.caps_directory, "masks", f"tpl-{template}")

    if args.roi_list is None:
        raise ValueError("A list of regions must be given.")
    else:
        check_mask_list(masks_location, args.roi_list, args.mask_pattern, cropping)

    for subject in os.listdir(path.join(args.caps_directory, "subjects")):
        subject_path = path.join(args.caps_directory, "subjects", subject)
        for session in os.listdir(subject_path):
            try:
                input_path = clinica_file_reader(
                    [subject], [session], args.caps_directory, FILE_TYPE
                )[0]
                output_path = path.join(
                    subject_path,
                    session,
                    "deeplearning_prepare_data",
                    "roi_based",
                    output_folder,
                )
                extract_roi(
                    input_path,
                    output_path,
                    masks_location,
                    args.mask_pattern,
                    cropping,
                    roi_list=args.roi_list,
                    uncrop_output=args.uncrop_output,
                )
            except ClinicaCAPSError:
                print("Subject %s session %s was not treated." % (subject, session))
