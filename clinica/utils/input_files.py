"""This module contains dictionaries used in inputs.py::clinica_{file|group}_reader().

These dictionaries describe files to grab.
"""

""" T1w """

# BIDS

T1W_NII = {"pattern": "sub-*_ses-*_t1w.nii*", "description": "T1w MRI"}

# T1-FreeSurfer

T1_FS_WM = {
    "pattern": "t1/freesurfer_cross_sectional/sub-*_ses-*/mri/wm.seg.mgz",
    "description": "segmentation of white matter (mri/wm.seg.mgz).",
    "needed_pipeline": "t1-freesurfer",
}

T1_FS_BRAIN = {
    "pattern": "t1/freesurfer_cross_sectional/sub-*_ses-*/mri/brain.mgz",
    "description": " extracted brain from T1w MRI (mri/brain.mgz).",
    "needed_pipeline": "t1-freesurfer",
}

T1_FS_ORIG_NU = {
    "pattern": "t1/freesurfer_cross_sectional/sub-*_ses-*/mri/orig_nu.mgz",
    "description": "intensity normalized volume generated after correction for"
    " non-uniformity in FreeSurfer (mri/orig_nu.mgz).",
    "needed_pipeline": "t1-freesurfer",
}

T1_FS_LONG_ORIG_NU = {
    "pattern": "t1/long-*/freesurfer_longitudinal/sub-*_ses-*.long.sub-*_*/mri/orig_nu.mgz",
    "description": "intensity normalized volume generated after correction for non-uniformity in FreeSurfer (orig_nu.mgz) in longitudinal",
    "needed_pipeline": "t1-freesurfer and t1-freesurfer longitudinal",
}

T1_FS_WM_SURF_R = {
    "pattern": "t1/freesurfer_cross_sectional/sub-*_ses-*/surf/rh.white",
    "description": "right white matter/gray matter border surface (rh.white).",
    "needed_pipeline": "t1-freesurfer",
}

T1_FS_LONG_SURF_R = {
    "pattern": "t1/long-*/freesurfer_longitudinal/sub-*_ses-*.long.sub-*_*/surf/rh.white",
    "description": "right white matter/gray matter border surface (rh.white) generated with t1-freesurfer-longitudinal.",
    "needed_pipeline": "t1-freesurfer and t1-freesurfer longitudinal",
}

T1_FS_LONG_SURF_L = {
    "pattern": "t1/long-*/freesurfer_longitudinal/sub-*_ses-*.long.sub-*_*/surf/lh.white",
    "description": "left white matter/gray matter border surface (lh.white) generated with t1-freesurfer-longitudinal.",
    "needed_pipeline": "t1-freesurfer and t1-freesurfer longitudinal",
}

T1_FS_WM_SURF_L = {
    "pattern": "t1/freesurfer_cross_sectional/sub-*_ses-*/surf/lh.white",
    "description": "left white matter/gray matter border surface (lh.white).",
    "needed_pipeline": "t1-freesurfer",
}

T1_FS_DESTRIEUX = {
    "pattern": "t1/freesurfer_cross_sectional/sub-*_ses-*/mri/aparc.a2009s+aseg.mgz",
    "description": "Destrieux-based segmentation (mri/aparc.a2009s+aseg.mgz).",
    "needed_pipeline": "t1-freesurfer",
}

T1_FS_DESTRIEUX_PARC_L = {
    "pattern": "t1/freesurfer_cross_sectional/sub-*_ses-*/label/lh.aparc.a2009s.annot",
    "description": "left hemisphere surface-based Destrieux parcellation (label/lh.aparc.a2009s.annot).",
    "needed_pipeline": "t1-freesurfer",
}

T1_FS_LONG_DESTRIEUX_PARC_L = {
    "pattern": "t1/long-*/freesurfer_longitudinal/sub-*_ses-*.long.sub-*_*/label/lh.aparc.a2009s.annot",
    "description": "left hemisphere surface-based Destrieux parcellation (label/lh.aparc.a2009s.annot) generated with t1-freesurfer-longitudinal.",
    "needed_pipeline": "t1-freesurfer and t1-freesurfer longitudinal",
}

T1_FS_LONG_DESTRIEUX_PARC_R = {
    "pattern": "t1/long-*/freesurfer_longitudinal/sub-*_ses-*.long.sub-*_*/label/rh.aparc.a2009s.annot",
    "description": "right hemisphere surface-based Destrieux parcellation (label/rh.aparc.a2009s.annot) generated with t1-freesurfer-longitudinal.",
    "needed_pipeline": "t1-freesurfer and t1-freesurfer longitudinal",
}

T1_FS_DESTRIEUX_PARC_R = {
    "pattern": "t1/freesurfer_cross_sectional/sub-*_ses-*/label/rh.aparc.a2009s.annot",
    "description": "right hemisphere surface-based Destrieux parcellation (label/rh.aparc.a2009s.annot).",
    "needed_pipeline": "t1-freesurfer",
}

T1_FS_DESIKAN = {
    "pattern": "t1/freesurfer_cross_sectional/sub-*_ses-*/mri/aparc+aseg.mgz",
    "description": "Desikan-based segmentation (mri/aparc.a2009s+aseg.mgz).",
    "needed_pipeline": "t1-freesurfer",
}

T1_FS_DESIKAN_PARC_L = {
    "pattern": "t1/freesurfer_cross_sectional/sub-*_ses-*/label/lh.aparc.annot",
    "description": "left hemisphere surface-based Desikan parcellation (label/lh.aparc.annot).",
    "needed_pipeline": "t1-freesurfer",
}

T1_FS_DESIKAN_PARC_R = {
    "pattern": "t1/freesurfer_cross_sectional/sub-*_ses-*/label/rh.aparc.annot",
    "description": "right hemisphere surface-based Desikan parcellation (label/rh.aparc.annot).",
    "needed_pipeline": "t1-freesurfer",
}

# T1-FreeSurfer-Template
T1_FS_T_DESTRIEUX = {
    "pattern": "freesurfer_unbiased_template/sub-*_long-*/mri/aparc.a2009s+aseg.mgz",
    "description": "Destrieux-based segmentation (mri/aparc.a2009s+aseg.mgz) from unbiased template.",
    "needed_pipeline": "t1-freesurfer-longitudinal or t1-freesurfer-template",
}

# T1-FreeSurfer-Longitudinal-Correction
T1_FS_LONG_DESIKAN_PARC_L = {
    "pattern": "t1/long-*/freesurfer_longitudinal/sub-*_ses-*.long.sub-*_*/label/lh.aparc.annot",
    "description": "left hemisphere surface-based Desikan parcellation (label/lh.aparc.annot) generated with t1-freesurfer-longitudinal.",
    "needed_pipeline": "t1-freesurfer and t1-freesurfer-longitudinal",
}

T1_FS_LONG_DESIKAN_PARC_R = {
    "pattern": "t1/long-*/freesurfer_longitudinal/sub-*_ses-*.long.sub-*_*/label/rh.aparc.annot",
    "description": "right hemisphere surface-based Desikan parcellation (label/rh.aparc.annot) generated with t1-freesurfer-longitudinal.",
    "needed_pipeline": "t1-freesurfer and t1-freesurfer-longitudinal",
}

T1W_LINEAR = {
    "pattern": "*space-MNI152NLin2009cSym_res-1x1x1_T1w.nii.gz",
    "description": "T1w image registered in MNI152NLin2009cSym space using t1-linear pipeline",
    "needed_pipeline": "t1-linear",
}

T1W_LINEAR_CROPPED = {
    "pattern": "*space-MNI152NLin2009cSym_desc-Crop_res-1x1x1_T1w.nii.gz",
    "description": "T1W Image registered using t1-linear and cropped "
    "(matrix size 169×208×179, 1 mm isotropic voxels)",
    "needed_pipeline": "t1-linear",
}

T1W_EXTENSIVE = {
    "pattern": "*space-Ixi549Space_desc-SkullStripped_T1w.nii.gz",
    "description": "T1w image skull-stripped registered in Ixi549Space space using clinicaDL preprocessing pipeline",
    "needed_pipeline": "t1-extensive",
}

T1W_TO_MNI_TRANSFORM = {
    "pattern": "*space-MNI152NLin2009cSym_res-1x1x1_affine.mat",
    "description": "Transformation matrix from T1W image to MNI space using t1-linear pipeline",
    "needed_pipeline": "t1-linear",
}

# T1-Volume


def t1_volume_native_tpm(tissue_number):
    import os

    from .spm import INDEX_TISSUE_MAP

    information = {
        "pattern": os.path.join(
            "t1",
            "spm",
            "segmentation",
            "native_space",
            f"*_*_T1w_segm-{INDEX_TISSUE_MAP[tissue_number]}_probability.nii*",
        ),
        "description": f"Tissue probability map {INDEX_TISSUE_MAP[tissue_number]} in native space",
        "needed_pipeline": "t1-volume-tissue-segmentation",
    }
    return information


def t1_volume_dartel_input_tissue(tissue_number):
    import os

    from .spm import INDEX_TISSUE_MAP

    information = {
        "pattern": os.path.join(
            "t1",
            "spm",
            "segmentation",
            "dartel_input",
            f"*_*_T1w_segm-{INDEX_TISSUE_MAP[tissue_number]}_dartelinput.nii*",
        ),
        "description": f"Dartel input for tissue probability map {INDEX_TISSUE_MAP[tissue_number]} from T1w MRI",
        "needed_pipeline": "t1-volume-tissue-segmentation",
    }
    return information


def t1_volume_native_tpm_in_mni(tissue_number, modulation):
    import os

    from .spm import INDEX_TISSUE_MAP

    if modulation:
        pattern_modulation = "on"
        description_modulation = "with"
    else:
        pattern_modulation = "off"
        description_modulation = "without"
    information = {
        "pattern": os.path.join(
            "t1",
            "spm",
            "segmentation",
            "normalized_space",
            f"*_*_T1w_segm-{INDEX_TISSUE_MAP[tissue_number]}_space-Ixi549Space_modulated-{pattern_modulation}_probability.nii*",
        ),
        "description": (
            f"Tissue probability map {INDEX_TISSUE_MAP[tissue_number]} based on "
            f"native MRI in MNI space (Ixi549) {description_modulation} modulation."
        ),
        "needed_pipeline": "t1-volume-tissue-segmentation",
    }
    return information


def t1_volume_template_tpm_in_mni(group_label, tissue_number, modulation):
    import os

    from .spm import INDEX_TISSUE_MAP

    if modulation:
        pattern_modulation = "on"
        description_modulation = "with"
    else:
        pattern_modulation = "off"
        description_modulation = "without"
    information = {
        "pattern": os.path.join(
            "t1",
            "spm",
            "dartel",
            f"group-{group_label}",
            f"*_T1w_segm-{INDEX_TISSUE_MAP[tissue_number]}_space-Ixi549Space_modulated-{pattern_modulation}_probability.nii*",
        ),
        "description": (
            f"Tissue probability map {INDEX_TISSUE_MAP[tissue_number]} based "
            f"on {group_label} template in MNI space (Ixi549) {description_modulation} modulation."
        ),
        "needed_pipeline": "t1-volume",
    }
    return information


def t1_volume_deformation_to_template(group_label):
    import os

    information = {
        "pattern": os.path.join(
            "t1",
            "spm",
            "dartel",
            f"group-{group_label}",
            f"sub-*_ses-*_T1w_target-{group_label}_transformation-forward_deformation.nii*",
        ),
        "description": f"Deformation from native space to group template {group_label} space.",
        "needed_pipeline": "t1-volume-create-dartel",
    }
    return information


def t1_volume_i_th_iteration_group_template(group_label, i):
    import os

    information = {
        "pattern": os.path.join(
            f"group-{group_label}",
            "t1",
            f"group-{group_label}_iteration-{i}_template.nii*",
        ),
        "description": f"Iteration #{i} of Dartel template {group_label}",
        "needed_pipeline": "t1-volume or t1-volume-create-dartel",
    }
    return information


def t1_volume_final_group_template(group_label):
    import os

    information = {
        "pattern": os.path.join(
            f"group-{group_label}", "t1", f"group-{group_label}_template.nii*"
        ),
        "description": f"T1w template file of group {group_label}",
        "needed_pipeline": "t1-volume or t1-volume-create-dartel",
    }
    return information


""" DWI """

# BIDS

DWI_NII = {"pattern": "dwi/sub-*_ses-*_dwi.nii*", "description": "DWI NIfTI"}

DWI_JSON = {"pattern": "dwi/sub-*_ses-*_dwi.json", "description": "DWI JSON file"}

DWI_BVAL = {"pattern": "dwi/sub-*_ses-*_dwi.bval", "description": "bval files"}

DWI_BVEC = {"pattern": "dwi/*_dwi.bvec", "description": "bvec files"}

FMAP_PHASEDIFF_JSON = {
    "pattern": "fmap/sub-*_ses-*_phasediff.json",
    "description": "phasediff JSON file",
}

FMAP_PHASEDIFF_NII = {
    "pattern": "fmap/sub-*_ses-*_phasediff.nii*",
    "description": "phasediff NIfTI volume",
}

FMAP_MAGNITUDE1_NII = {
    "pattern": "fmap/sub-*_ses-*_magnitude1.nii*",
    "description": "magnitude1 file",
}

# CAPS

DWI_PREPROC_NII = {
    "pattern": "dwi/preprocessing/sub-*_ses-*_dwi_space-*_preproc.nii*",
    "description": "preprocessed DWI",
    "needed_pipeline": "dwi-preprocessing-using-t1 or dwi-preprocessing-using-fieldmap",
}

DWI_PREPROC_BRAINMASK = {
    "pattern": "dwi/preprocessing/sub-*_ses-*_dwi_space-*_brainmask.nii*",
    "description": "b0 brainmask",
    "needed_pipeline": "dwi-preprocessing-using-t1 or dwi-preprocessing-using-fieldmap",
}

DWI_PREPROC_BVEC = {
    "pattern": "dwi/preprocessing/sub-*_ses-*_dwi_space-*_preproc.bvec",
    "description": "preprocessed bvec",
    "needed_pipeline": "dwi-preprocessing-using-t1 or dwi-preprocessing-using-fieldmap",
}

DWI_PREPROC_BVAL = {
    "pattern": "dwi/preprocessing/*_dwi_space-*_preproc.bval",
    "description": "preprocessed bval",
    "needed_pipeline": "dwi-preprocessing-using-t1 or dwi-preprocessing-using-fieldmap",
}

""" PET """

# BIDS


def bids_pet_nii(acq_label):
    import os

    information = {
        "pattern": os.path.join("pet", f"sub-*_ses-*_acq-{acq_label}_pet.nii*"),
        "description": f"PET data with {acq_label} tracer",
    }
    return information


# PET-Volume


def pet_volume_normalized_suvr_pet(
    acq_label,
    group_label,
    suvr_reference_region,
    use_brainmasked_image,
    use_pvc_data,
    fwhm=0,
):
    import os

    if use_brainmasked_image:
        mask_key_value = "_mask-brain"
        mask_description = "brain-masked"
    else:
        mask_key_value = ""
        mask_description = "full"

    if use_pvc_data:
        pvc_key_value = "_pvc-rbv"
        pvc_description = "using RBV method for PVC"
    else:
        pvc_key_value = ""
        pvc_description = "without PVC"

    if fwhm:
        fwhm_key_value = f"_fwhm-{fwhm}mm"
        fwhm_description = f"with {fwhm}mm smoothing"
    else:
        fwhm_key_value = ""
        fwhm_description = "with no smoothing"

    suvr_key_value = f"_suvr-{suvr_reference_region}"

    information = {
        "pattern": os.path.join(
            "pet",
            "preprocessing",
            f"group-{group_label}",
            f"*_acq-{acq_label}_pet"
            f"_space-Ixi549Space{pvc_key_value}{suvr_key_value}{mask_key_value}{fwhm_key_value}_pet.nii*",
        ),
        "description": (
            f"{mask_description} SUVR map (using {suvr_reference_region} region) of {acq_label}-PET "
            f"{pvc_description} and {fwhm_description} in Ixi549Space space based on {group_label} DARTEL template"
        ),
        "needed_pipeline": "pet-volume",
    }
    return information


# pet-linear


def pet_linear_nii(acq_label, suvr_reference_region, uncropped_image):
    import os

    if uncropped_image:
        description = ""
    else:
        description = "_desc-Crop"

    information = {
        "pattern": os.path.join(
            "pet_linear",
            f"*_acq-{acq_label}_pet_space-MNI152NLin2009cSym{description}_res-1x1x1_suvr-{suvr_reference_region}_pet.nii.gz",
        ),
        "description": "",
        "needed_pipeline": "pet-linear",
    }
    return information
