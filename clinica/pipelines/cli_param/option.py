"""Common CLI options used by Clinica pipelines."""

import click

from clinica.utils.pet import LIST_SUVR_REFERENCE_REGIONS

from .option_group import option

acq_label = option(
    "-acq",
    "--acq_label",
    help="Name of the label given to the PET acquisition, specifying the tracer used (acq-<acq_label>).",
)

dartel_tissues = option(
    "-dt",
    "--dartel_tissues",
    multiple=True,
    type=click.IntRange(1, 6),
    default=(1, 2, 3),
    show_default=True,
    help="Tissues to use for DARTEL template calculation.",
)

dont_save_warped_unmodulated = option(
    "-dswu",
    "--dont_save_warped_unmodulated",
    is_flag=True,
    help="Do not save warped unmodulated images for tissues specified in --tissue_classes flag.",
)

initrand = option(
    "--initrand",
    default=False,
    is_flag=True,
    show_default=True,
    help="Initialize the random number generator used when estimating hyperparameters in FSL eddy.",
)

low_bval = option(
    "--low_bval",
    default=5,
    show_default=True,
    help="Define the b0 volumes as all volume bval <= low_bval.",
)

modulate = option(
    "--modulate/--no_modulate",
    "modulate",
    default=True,
    show_default=True,
    help="Modulate output images, no modulation preserves concentrations.",
)

n_procs = option(
    "-np",
    "--n_procs",
    type=int,
    help="Number of cores used to run in parallel.",
)

overwrite_outputs = option(
    "-overwrite",
    "--overwrite_outputs",
    is_flag=True,
    help="Force overwrite of output files in CAPS folder.",
)

pvc_psf_tsv = option(
    "-psf",
    "--pvc_psf_tsv",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
    help=(
        "TSV file containing for each PET image its point spread function (PSF) measured "
        "in mm at x, y & z coordinates. Columns must contain: "
        "participant_id, session_id, acq_label, psf_x, psf_y and psf_z."
    ),
)

save_warped_modulated = option(
    "-swm",
    "--save_warped_modulated",
    is_flag=True,
    help="Save warped modulated images for tissues specified in --tissue_classes flag.",
)

smooth = option(
    "-s",
    "--smooth",
    multiple=True,
    default=[8],
    show_default=True,
    help="Specify the different isomorphic FWHM in millimeters to smooth the image.",
)

subjects_sessions_tsv = option(
    "-tsv",
    "--subjects_sessions_tsv",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
    help="TSV file containing a list of subjects with their sessions.",
)

suvr_reference_region = option(
    "-suvr",
    "--suvr_reference_region",
    type=click.Choice(LIST_SUVR_REFERENCE_REGIONS),
    help=(
        "Intensity normalization using the average PET uptake in reference regions "
        "resulting in a standardized uptake value ratio (SUVR) map. It can be "
        "cerebellumPons (used for amyloid tracers) or pons (used for 18F-FDG tracers)."
    ),
)

tissues = option(
    "-t",
    "--tissues",
    multiple=True,
    type=click.IntRange(1, 6),
    default=(1, 2, 3),
    show_default=True,
    help="Tissues to create flow fields to DARTEL template.",
)

tissue_classes = option(
    "-tc",
    "--tissue_classes",
    multiple=True,
    type=click.IntRange(1, 6),
    default=(1, 2, 3),
    show_default=True,
    help=(
        "Tissue classes (1: gray matter (GM), 2: white matter (WM), 3: cerebrospinal fluid (CSF), 4: bone, "
        "5: soft-tissue, 6: background) to save."
    ),
)

tissue_probability_maps = option(
    "-tpm",
    "--tissue_probability_maps",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
    help="Tissue probability maps to use for segmentation.",
)

use_cuda = option(
    "--use_cuda",
    default=False,
    is_flag=True,
    show_default=True,
    help="Use the CUDA implementation of FSL eddy. See https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy/UsersGuide#The_eddy_executables for details.",
)

use_pvc_data = option(
    "-pvc",
    "--use_pvc_data",
    is_flag=True,
    help="Use PET data with partial value correction.",
)

voxel_size = option(
    "-vs",
    "--voxel_size",
    type=click.Tuple([float, float, float]),
    default=(1.5, 1.5, 1.5),
    show_default=True,
    help="Specifying the voxel sizes of the output image.",
)

working_directory = option(
    "-wd",
    "--working_directory",
    type=click.Path(file_okay=False, writable=True, resolve_path=True),
    help="Temporary directory to store pipelines intermediate results.",
)

yes = option(
    "-y",
    "--yes",
    is_flag=True,
    help="Execute the pipeline even if input images are not centered without asking for more user input.",
)

atlas_path = option(
    "-ap",
    "--atlas_path",
    type=click.Path(exists=True, dir_okay=True, resolve_path=True),
    help="Compute atlases at the end of the path",
)

delete_cache = option(
    "-dc",
    "--delete_cache",
    is_flag=True,
    help="If True, large intermediary files will be deleted before the end of the pipeline.",  # This option is only implemented in dwi_preprocessing_using_t1 for now.
)
