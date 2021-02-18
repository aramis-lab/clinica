# coding: utf8


def extract_slices(input_tensor, slice_direction=0, slice_mode="single"):
    """Extracts the slices from three directions

    This function extracts slices form the preprocesed nifti image.  The
    direction of extraction can be defined either on sagital direction (0),
    cornal direction (1) or axial direction (other). The output slices can be
    stores following two modes: single (1 channel) ou RGB (3 channels, all the
    same).


    Args:
        input_tensor: tensor version of the nifti MRI.
        slice_direction: which axis direction that the slices were extracted
        slice_mode: 'single' or 'RGB'.

    Returns:
        file: multiple tensors saved on the disk, suffixes corresponds to
            indexes of the slices. Same location than input file.
    """
    import os

    import torch

    image_tensor = torch.load(input_tensor)
    # reshape the tensor, delete the first dimension for slice-level
    image_tensor = image_tensor.view(
        image_tensor.shape[1], image_tensor.shape[2], image_tensor.shape[3]
    )

    # sagital
    # M and N correspond to the first and last slices (if need to remove)
    M = 0
    N = 0
    slice_list_sag = range(
        M, image_tensor.shape[0] - N
    )  # delete the first M slices and last N slices

    basedir = os.getcwd()
    input_tensor_filename = os.path.basename(input_tensor)

    txt_idx = input_tensor_filename.rfind("_")
    it_filename_prefix = input_tensor_filename[0:txt_idx]
    it_filename_suffix = input_tensor_filename[txt_idx:]

    output_file_original = []
    output_file_rgb = []
    if slice_direction == 0:
        for index_slice, index_slice_list in zip(
            slice_list_sag, range(len(slice_list_sag))
        ):
            # for i in slice_list:
            # sagital
            slice_select_sag = image_tensor[index_slice, :, :]

            extracted_slice_original_sag = slice_select_sag.unsqueeze(
                0
            )  # shape should be 1 * W * L

            # train for transfer learning, creating the fake RGB image.
            slice_select_sag = (slice_select_sag - slice_select_sag.min()) / (
                slice_select_sag.max() - slice_select_sag.min()
            )
            extracted_slice_rgb_sag = torch.stack(
                (slice_select_sag, slice_select_sag, slice_select_sag)
            )  # shape should be 3 * W * L

            # save into .pt format
            if slice_mode == "single":
                output_file_original.append(
                    os.path.join(
                        basedir,
                        it_filename_prefix
                        + "_axis-sag_channel-single_slice-"
                        + str(index_slice)
                        + it_filename_suffix,
                    )
                )
                torch.save(
                    extracted_slice_original_sag.clone(),
                    output_file_original[index_slice_list],
                )
            elif slice_mode == "rgb":
                output_file_rgb.append(
                    os.path.join(
                        basedir,
                        it_filename_prefix
                        + "_axis-sag_channel-rgb_slice-"
                        + str(index_slice)
                        + it_filename_suffix,
                    )
                )
                torch.save(
                    extracted_slice_rgb_sag.clone(), output_file_rgb[index_slice_list]
                )

    elif slice_direction == 1:
        # cornal
        slice_list_cor = range(
            M, image_tensor.shape[1] - N
        )  # delete the first M slices and last N slices
        for index_slice, index_slice_list in zip(
            slice_list_cor, range(len(slice_list_cor))
        ):
            # for i in slice_list:
            # sagital
            slice_select_cor = image_tensor[:, index_slice, :]

            extracted_slice_original_cor = slice_select_cor.unsqueeze(
                0
            )  # shape should be 1 * W * L

            # train for transfer learning, creating the fake RGB image.
            slice_select_cor = (slice_select_cor - slice_select_cor.min()) / (
                slice_select_cor.max() - slice_select_cor.min()
            )
            extracted_slice_rgb_cor = torch.stack(
                (slice_select_cor, slice_select_cor, slice_select_cor)
            )  # shape should be 3 * W * L

            # save into .pt format
            if slice_mode == "single":
                output_file_original.append(
                    os.path.join(
                        basedir,
                        it_filename_prefix
                        + "_axis-cor_channel-single_slice-"
                        + str(index_slice)
                        + it_filename_suffix,
                    )
                )
                torch.save(
                    extracted_slice_original_cor.clone(),
                    output_file_original[index_slice_list],
                )
            elif slice_mode == "rgb":
                output_file_rgb.append(
                    os.path.join(
                        basedir,
                        it_filename_prefix
                        + "_axis-cor_channel-rgb_slice-"
                        + str(index_slice)
                        + it_filename_suffix,
                    )
                )
                torch.save(
                    extracted_slice_rgb_cor.clone(), output_file_rgb[index_slice_list]
                )

    else:

        # axial
        slice_list_axi = range(
            M, image_tensor.shape[2] - N
        )  # delete the first M slices and last N slices
        for index_slice, index_slice_list in zip(
            slice_list_axi, range(len(slice_list_axi))
        ):
            # for i in slice_list:
            # sagital
            slice_select_axi = image_tensor[:, :, index_slice]

            extracted_slice_original_axi = slice_select_axi.unsqueeze(
                0
            )  # shape should be 1 * W * L

            # train for transfer learning, creating the fake RGB image.
            slice_select_axi = (slice_select_axi - slice_select_axi.min()) / (
                slice_select_axi.max() - slice_select_axi.min()
            )
            extracted_slice_rgb_axi = torch.stack(
                (slice_select_axi, slice_select_axi, slice_select_axi)
            )  # shape should be 3 * W * L

            # save into .pt format
            if slice_mode == "single":
                output_file_original.append(
                    os.path.join(
                        basedir,
                        it_filename_prefix
                        + "_axis-axi_channel-single_slice-"
                        + str(index_slice)
                        + it_filename_suffix,
                    )
                )
                torch.save(
                    extracted_slice_original_axi.clone(),
                    output_file_original[index_slice_list],
                )
            elif slice_mode == "rgb":
                output_file_rgb.append(
                    os.path.join(
                        basedir,
                        it_filename_prefix
                        + "_axis-axi_channel-rgb_slice-"
                        + str(index_slice)
                        + it_filename_suffix,
                    )
                )
                torch.save(
                    extracted_slice_rgb_axi.clone(), output_file_rgb[index_slice_list]
                )

    return output_file_rgb, output_file_original


def extract_patches(input_tensor, patch_size, stride_size):
    """Extracts the patches

    This function extracts patches form the preprocesed nifti image. Patch size
    if provieded as input and also the stride size. If stride size is smaller
    than the patch size an overlap exist between consecutive patches. If stride
    size is equal to path size there is no overlap. Otherwise, unprocessed
    zones can exits.

    Args:
        input_tensor: tensor version of the nifti MRI.
        patch_size: size of a single patch.
        stride_size: size of the stride leading to next patch.

    Returns:
        file: multiple tensors saved on the disk, suffixes corresponds to
            indexes of the patches. Same location than input file.
    """
    import os

    import torch

    basedir = os.getcwd()
    image_tensor = torch.load(input_tensor)

    # use classifiers tensor.upfold to crop the patch.
    patches_tensor = (
        image_tensor.unfold(1, patch_size, stride_size)
        .unfold(2, patch_size, stride_size)
        .unfold(3, patch_size, stride_size)
        .contiguous()
    )
    # the dimension of patch_tensor should be [1, patch_num1, patch_num2, patch_num3, patch_size1, patch_size2, patch_size3]
    patches_tensor = patches_tensor.view(-1, patch_size, patch_size, patch_size)

    input_tensor_filename = os.path.basename(input_tensor)
    txt_idx = input_tensor_filename.rfind("_")
    it_filename_prefix = input_tensor_filename[0:txt_idx]
    it_filename_suffix = input_tensor_filename[txt_idx:]

    output_patch = []
    for index_patch in range(patches_tensor.shape[0]):
        extracted_patch = patches_tensor[index_patch, ...].unsqueeze_(
            0
        )  # add one dimension
        # save into .pt format
        output_patch.append(
            os.path.join(
                basedir,
                it_filename_prefix
                + "_patchsize-"
                + str(patch_size)
                + "_stride-"
                + str(stride_size)
                + "_patch-"
                + str(index_patch)
                + it_filename_suffix,
            )
        )
        torch.save(extracted_patch.clone(), output_patch[index_patch])

    return output_patch


def save_as_pt(input_img):
    """Saves PyTorch tensor version of the nifti image

    This function convert nifti image to tensor (.pt) version of the image.
    Tensor version is saved at the same location than input_img.

    Args:
        input_tensor: tensor version of the nifti MRI.

    Returns:
        filename (str): single tensor file  saved on the disk. Same location than input file.

    """

    import os

    import nibabel as nib
    import torch

    basedir = os.getcwd()
    image_array = nib.load(input_img).get_fdata()
    image_tensor = torch.from_numpy(image_array).unsqueeze(0).float()
    # make sure the tensor dtype is torch.float32
    output_file = os.path.join(
        basedir, os.path.basename(input_img).split(".nii.gz")[0] + ".pt"
    )
    # save
    torch.save(image_tensor.clone(), output_file)

    return output_file
