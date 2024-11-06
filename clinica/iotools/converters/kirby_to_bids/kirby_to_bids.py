"""Convert KIRBY dataset (https://brain-development.org/ixi-dataset/) to BIDS."""

from pathlib import Path
from typing import Optional

import nibabel as nb
import numpy as np
import csv

from clinica.iotools.bids_utils import write_modality_agnostic_files
from clinica.iotools.converters.kirby_to_bids.kirby_to_bids_utils import (
  create_bids_structure,
  convert_to_nii_gz,
  find_nii_files,
  normalize_dashes,
  replace_dashes_with_underscore,
)
from clinica.utils.filemanip import UserProvidedPath

__all__ = ["convert"]

import os
import pandas as pd

def convert(
    path_to_dataset: str,
    bids_dir: str,
    path_to_clinical: str,
    subjects: Optional[str] = None,
    n_procs: Optional[int] = 1,
    **kwargs,
):
    """Convert KIRBY data to BIDS format from Excel clinical data."""
    clinical_data_file = None
    for file in os.listdir(path_to_clinical):
        if file.endswith('.xlsx'):
            clinical_data_file = os.path.join(path_to_clinical, file)
            break

    if not clinical_data_file:
        raise FileNotFoundError(f"No clinical data Excel file found in {path_to_clinical}")

    clinical_data = pd.read_excel(clinical_data_file)
    clinical_data_filtered = clinical_data[['MPRAGE', 'Age', 'Sex', 'Fiducial', 'Subject ID', 'Visit ID']]
    clinical_data_filtered['MPRAGE'] = clinical_data_filtered['MPRAGE'].apply(replace_dashes_with_underscore).str.strip()
    clinical_data_filtered.reset_index(drop=True, inplace=True)

    participants_data = {}
    session_count = {}
    
    # Get all .nii files (directly or within subfolders)
    nii_files = find_nii_files(path_to_dataset)

    # Traverse found nii files
    for file_path in nii_files:
        file = os.path.basename(file_path).replace('.nii', '').strip()  # Remove file extension and strip whitespaces
        file_normalized = replace_dashes_with_underscore(file)
        
        if 'MPRAGE' in file_normalized:
            # Find the exact match in the clinical data's MPRAGE column
            clinical_row = clinical_data_filtered[clinical_data_filtered['MPRAGE'] == file_normalized]

            # Check if any matching rows were found
            if clinical_row.empty:
                print(f"No matching clinical data found for file: {file_normalized}")
                continue

            # Extract the first matching row (in case multiple matches are found)
            clinical_row = clinical_row.iloc[0]

            # Extract relevant clinical information
            subject_id = clinical_row['Subject ID']
            session_id = clinical_row['Visit ID']
            age = clinical_row['Age']
            sex = clinical_row['Sex']
            handedness = clinical_row['Fiducial']

            # Keep track of how many sessions have been processed for this subject
            if subject_id not in session_count:
                session_count[subject_id] = 1
            else:
                session_count[subject_id] += 1

            # Assign session as "ses-01", "ses-02", etc.
            session_label = f"{session_count[subject_id]:02d}"

            # Create BIDS structure and move the file
            create_bids_structure(subject_id, session_label, file_path, bids_dir)

            # Store participant data (with baseline age if needed)
            if subject_id not in participants_data:
                participants_data[subject_id] = {'participant_id': f"sub-KKI{subject_id}",
                                                 'sex': sex,
                                                 'age': age,
                                                 'handedness': handedness}

    # Write participants.csv
    participants_csv = os.path.join(bids_dir, 'participants.csv')
    with open(participants_csv, 'w', newline='') as participants_file:
        participants_writer = csv.writer(participants_file)
        participants_writer.writerow(['participant_id', 'sex', 'age', 'handedness'])

        for participant_info in participants_data.values():
            participants_writer.writerow([participant_info['participant_id'],
                                          participant_info['sex'],
                                          participant_info['age'],
                                          participant_info['handedness']])

    # Write sessions.tsv for each subject
    for subject_id, session_count in session_count.items():
        sessions_file = os.path.join(bids_dir, f"sub-KKI{subject_id}", 'sessions.tsv')
        os.makedirs(os.path.dirname(sessions_file), exist_ok=True)

        with open(sessions_file, 'w', newline='') as session_file:
            session_writer = csv.writer(session_file, delimiter='\t')
            session_writer.writerow(['session_id', 'age'])

            for session_num in range(1, session_count + 1):
                session_writer.writerow([f"ses-{session_num:02d}", participants_data[subject_id]['age']])

    print(f"KIRBY BIDS conversion completed using the clinical data from {clinical_data_file}.")
