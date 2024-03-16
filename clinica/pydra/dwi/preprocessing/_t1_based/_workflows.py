from pydra import Workflow

__all__ = ["build_b0_flirt_pipeline"]


def build_b0_flirt_pipeline(num_b0s: int, name: str = "b0_coregistration") -> Workflow:
    """Rigid registration of the B0 dataset onto the first volume.

    Rigid registration is achieved using FLIRT and the normalized correlation.

    Parameters
    ----------
    num_b0s : int
        Number of B0 volumes in the dataset.

    name : str, optional
        Name of the workflow.
        Default="b0_coregistration".

    Inputnode:
        in_file(str): B0 dataset.

    Outputnode
        out_b0_reg(str): The set of B0 volumes registered to the first volume.
    """
    from pydra.tasks import fsl

    from ._tasks import merge_nifti_images_in_time_dimension_task

    wf = Workflow(name, input_spec=["dwi_image_small_b_volumes"])
    wf.add(
        fsl.ROI(
            name="b0_reference",
            input_image=wf.lzin.dwi_image_small_b_volumes,
            tmin=0,
            tsize=1,
        )
    )
    wf.add(
        fsl.ROI(
            name="b0_moving",
            input_image=wf.lzin.dwi_image_small_b_volumes,
            tmin=1,
            tsize=num_b0s - 1,
        )
    )
    wf.add(
        fsl.Split(
            name="split_b0_moving",
            input_image=wf.b0_moving.lzout.output_image,
            direction="t",
        )
    )
    wf.add(
        fsl.BET(
            name="bet_ref",
            input_image=wf.b0_reference.lzout.output_image,
            fractional_intensity_threshold=0.3,
            save_brain_mask=True,
            with_robust_brain_center_estimation=True,
        )
    )
    wf.add(
        fsl.maths.FSLMaths(
            name="mask_dilate",
            input_image=wf.bet_ref.lzout.output_image,
            nan2zeros=True,  # TODO: fix this
            args="-kernel sphere 5 -dilM",  # TODO: fix this
        )
    )
    wf.add(
        fsl.FLIRT(
            name="b0_co_registration",
            input_image=wf.split_b0_moving.lzout.output_image,
            reference_image=wf.b0_reference.lzout.output_image,
            reference_weights=wf.mask_dilate.lzout.output_image,
            input_weights=wf.mask_dilate.lzout.output_image,
            interpolation="spline",
            degrees_of_freedom=6,
            num_bins=50,
            verbose=True,
            cost_function="corratio",
            padding_size=10,  # TODO: fix
            search_range_x=[-4, 4],
            search_range_y=[-4, 4],
            search_range_z=[-4, 4],
            fine_search=1,  # TODO: fix
            coarse_search=10,  # TODO: fix
        )
    ).split("input_image")
    wf.add(
        fsl.maths.Threshold(
            name="remove_negative",
            input_image=wf.b0_co_registration.lzout.out_image,
            threshold=0.0,
        )
    )
    wf.add(
        fsl.Merge(
            name="merge_registered_b0s",
            input_images=wf.remove_negative.lzout.output_image,
            dimension="t",
        )
    )
    wf.add(
        merge_nifti_images_in_time_dimension_task(
            image1=wf.b0_reference.lzout.output_image,
            image2=wf.merge_registered_b0s.lzout.output_image,
        )
    )
    wf.set_output(
        {
            "out_file": wf.merge_nifti_images_in_time_dimension_task.lzout.output_image,
            "out_xfms": wf.b0_co_registration.lzout.output_matrix,
        }
    )

    return wf
