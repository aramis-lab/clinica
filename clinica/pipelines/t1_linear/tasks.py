def run_n4biasfieldcorrection_task(
    input_image: str,
    bspline_fitting_distance: int,
    output_prefix=None,
    output_dir=None,
    save_bias=False,
    verbose=False,
) -> str:
    from pathlib import Path

    from clinica.pipelines.t1_linear.anat_linear_utils import run_n4biasfieldcorrection

    if output_dir:
        output_dir = Path(output_dir)

    return str(
        run_n4biasfieldcorrection(
            Path(input_image),
            bspline_fitting_distance,
            output_prefix,
            output_dir,
            save_bias,
            verbose,
        )
    )


def run_ants_registration_task(
    fixed_image: str,
    moving_image: str,
    random_seed: int,
    output_prefix=None,
    output_dir=None,
) -> tuple:
    from pathlib import Path

    from clinica.pipelines.t1_linear.anat_linear_utils import run_ants_registration

    if output_dir:
        output_dir = Path(output_dir)

    warped_image_output_path, transformation_matrix_output_path = run_ants_registration(
        Path(fixed_image),
        Path(moving_image),
        random_seed,
        output_prefix,
        output_dir,
    )

    return str(warped_image_output_path), str(transformation_matrix_output_path)
