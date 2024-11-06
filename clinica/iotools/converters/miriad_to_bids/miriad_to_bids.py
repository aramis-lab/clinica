"""Convert MIRIAD dataset to BIDS."""

from pathlib import Path
from typing import Optional

import os
import shutil
import csv
import pandas as pd
from clinica.utils.filemanip import UserProvidedPath

def convert(
    path_to_dataset: str,
    bids_dir: str,
    path_to_clinical: str,
    subjects: Optional[str] = None,
    n_procs: Optional[int] = 1,
    **kwargs,
):
    """Convert MIRIAD data to BIDS format without removing original .nii files."""
    from clinica.iotools.converters.miriad_to_bids.miriad_to_bids_utils import create_bids_structure, parse_filename, convert_to_nii_gz
    metadata_csv = 'metadata.csv'
    
    # Load clinical data
    clinical_data_file = None
    for file in os.listdir(path_to_clinical):
        if file.endswith('.csv'):
            clinical_data_file = os.path.join(path_to_clinical, file)
            break

    if not clinical_data_file:
        raise FileNotFoundError(f"No clinical data CSV found in {path_to_clinical}")

    clinical_data = pd.read_csv(clinical_data_file)

    # Prepare CSV
    with open(metadata_csv, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['cohort', 'subject_id', 'diagnosis', 'gender', 'session', 'run', 'input_file', 'output_file'])
        
        participants_data = {}
        sessions_data = []

        # Traverse the input directory
        for root, dirs, files in os.walk(path_to_dataset):
            for file in files:
                if file.endswith('.nii'):
                    # Extract information from filename
                    parts = file.split('_')
                    cohort = parts[0]   # miriad
                    subject_id = parts[1]  # 215
                    diagnosis = parts[2]  # AD (Alzheimer's) or HC (Healthy Control)
                    gender = parts[3]     # M or F
                    session = parts[4].lstrip('0')  # Session number
                    run_number = parts[6].replace('.nii', '')  # Scan number from MR_1 or MR_2

                    bids_subject_id = f"sub-{subject_id}"
                    bids_session_id = f"ses-{session}"
                    
                    # Original file path
                    original_file_path = os.path.join(root, file)
                    
                    # Extract MR ID
                    mr_id = f"{cohort}_{subject_id}_{session}_MR_{run_number}"

                    # Extract relevant clinical information from the clinical data
                    clinical_row = clinical_data[clinical_data['MR ID'] == mr_id]
                    if clinical_row.empty:
                        print(f"Clinical data not found for MR ID: {mr_id}")
                        continue

                    age = clinical_row['Age'].values[0]
                    group = clinical_row['Group'].values[0]  # HC or AD
                    gender_clinical = clinical_row['M/F'].values[0]  # M or F

                    # Write metadata CSV
                    csvwriter.writerow([cohort, subject_id, diagnosis, gender, session, run_number, original_file_path, bids_subject_id])

                    # Track baseline age (minimum age for each subject)
                    if subject_id not in participants_data or participants_data[subject_id]['age'] > age:
                        participants_data[subject_id] = {
                            'participant_id': f"sub-MIRIAD{subject_id}", 
                            'sex': gender_clinical, 
                            'diagnosis': group, 
                            'age': age
                        }

                    # Prepare sessions data
                    sessions_data.append([f"sub-MIRIAD{subject_id}", f"ses-{session}", age])

                    # Create BIDS structure and copy file with run number
                    create_bids_structure(subject_id, session, run_number, cohort, diagnosis, gender, original_file_path, path_to_dataset, bids_dir, path_to_clinical)

    # Write participants.csv with baseline age (minimum age for each subject)
    participants_csv = os.path.join(bids_dir, 'participants.csv')
    with open(participants_csv, 'w', newline='') as participants_file:
        participants_writer = csv.writer(participants_file)
        participants_writer.writerow(['participant_id', 'sex', 'diagnosis', 'age'])
        
        # Write the baseline age (minimum age) for each subject
        for participant_info in participants_data.values():
            participants_writer.writerow([participant_info['participant_id'], 
                                        participant_info['sex'], 
                                        participant_info['diagnosis'], 
                                        participant_info['age']])

    # Write sessions.tsv for each subject
    subject_sessions = {}
    for session in sessions_data:
        subject_id, session_id, age = session
        if subject_id not in subject_sessions:
            subject_sessions[subject_id] = []
        subject_sessions[subject_id].append([session_id, age])

    for subject_id, sessions in subject_sessions.items():
        sessions_file = os.path.join(bids_dir, subject_id, 'sessions.tsv')
        os.makedirs(os.path.dirname(sessions_file), exist_ok=True)
        
        with open(sessions_file, 'w', newline='') as session_file:
            session_writer = csv.writer(session_file, delimiter='\t')
            session_writer.writerow(['session_id', 'age'])
            session_writer.writerows(sessions)

    print(f"BIDS conversion completed, clinical data loaded from {clinical_data_file}.")
