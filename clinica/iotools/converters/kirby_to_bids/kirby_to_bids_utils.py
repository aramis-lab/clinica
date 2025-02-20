import json
import re
import shutil
from enum import Enum
from pathlib import Path
from typing import List, Optional, Union
import os

import nibabel as nib
import pandas as pd
from nilearn.image import concat_imgs

from clinica.iotools.bids_utils import StudyName, bids_id_factory
from clinica.utils.stream import cprint, log_and_raise

__all__ = [
    "read_clinical_data",
    "define_participants",
    "write_subject_data",
    "write_sessions",
    "write_scans",
    "write_participants",
    "check_modalities",
]


# Function to recursively find all files with .nii extension in input directory
def find_nii_files(directory):
    nii_files = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.nii') and 'MPRAGE' in file:
                nii_files.append(os.path.join(root, file))
    return nii_files

# Function to normalize dashes (replace any type of dash with a standard hyphen '-')
def normalize_dashes(text):
    if isinstance(text, str):
        return re.sub(r'[\u2013\u2014\u2212]', '-', text)  # Replaces en dash, em dash, and other similar symbols
    return text

def convert_to_nii_gz(input_file):
    """Convert a .nii file to .nii.gz format without deleting the original .nii file."""
    if input_file.endswith(".nii.gz"):
        return input_file
    img = nib.load(input_file)
    output_file = input_file.replace(".nii", ".nii.gz")
    nib.save(img, output_file)
    return output_file

def replace_dashes_with_underscore(string):
    """Replace dashes with underscores for standardized comparison."""
    return string.replace("-", "_")

# Helper function to create BIDS folders and move files
def create_bids_structure(subject_id, session, input_file, output_dir):
    """Create BIDS folder structure and move files into it."""
    sub_id = f"sub-KKI{subject_id}"
    ses_id = f"ses-{session}"

    # Create output directory for this subject/session
    anat_dir = os.path.join(output_dir, sub_id, ses_id, 'anat')
    os.makedirs(anat_dir, exist_ok=True)
    
    # Convert the input file to .nii.gz if necessary
    input_file_gz = convert_to_nii_gz(input_file)
    
    # Destination filename in BIDS format
    bids_filename = f"{sub_id}_{ses_id}_T1w.nii.gz"
    
    # Copy and rename the file to BIDS format
    shutil.copy(input_file_gz, os.path.join(anat_dir, bids_filename))