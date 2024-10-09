"""Convert MIRIAD dataset to BIDS."""

from pathlib import Path
from typing import Optional

import os
import shutil
import csv
import pandas as pd
from clinica.utils.filemanip import UserProvidedPath

def convert(
    path_to_dataset: UserProvidedPath,
    bids_dir: UserProvidedPath,
    path_to_clinical: UserProvidedPath,
    subjects: Optional[UserProvidedPath] = None,
    n_procs: Optional[int] = 1,
    **kwargs,
):
    """_summary_

    Args:
        path_to_dataset (UserProvidedPath): _description_
        bids_dir (UserProvidedPath): _description_
        path_to_clinical (UserProvidedPath): _description_
        subjects (Optional[UserProvidedPath], optional): _description_. Defaults to None.
        n_procs (Optional[int], optional): _description_. Defaults to 1.
    """
    from clinica.iotools.converters.miriad_to_bids.miriad_to_bids_utils import create_bids_structure, parse_filename
    metadata_csv = 'metadata.csv'
    # Load clinical data
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
        csvwriter.writerow(['cohort', 'subject_id', 'diagnosis', 'gender', 'session', 'input_file', 'output_file'])
        
        participants_data = {}
        sessions_data = []

        # Traverse the input directory
        for root, dirs, files in os.walk(path_to_dataset):
            for file in files:
                if file.endswith('.nii'):
                    # Example: miriad_215_AD_M_01_MR_1.nii
                    parts = file.split('_')
                    
                    # Extract information from filename
                    cohort = parts[0]   # miriad
                    subject_id = parts[1]  # 215
                    diagnosis = parts[2]  # AD (Alzheimer's) or HC (Healthy Control)
                    gender = parts[3]     # M or F
                    session = parts[4].lstrip('0')    # Session number
                    session_alt = parts[4].lstrip('0')    # Session number
                    scan_number = parts[6].replace('.nii', '')  # Scan number from MR_1 or MR_2

                        # Parse subject ID, session ID, and run ID from the filename
           # subject_id, session_id, run_id = parse_filename(file)

                    # Extract MR ID
                    mr_id = f"{cohort}_{subject_id}_{session}_MR_{scan_number}"

                    # Extract relevant clinical information from the clinical data
                    clinical_row = clinical_data[clinical_data['MR ID'] == mr_id]
                    if clinical_row.empty:
                        print(f"Clinical data not found for MR ID: {mr_id}")
                        continue

                    age = clinical_row['Age'].values[0]
                    group = clinical_row['Group'].values[0]  # HC or AD
                    gender_clinical = clinical_row['M/F'].values[0]  # M or F

                    # Full path of input file
                    input_file = os.path.join(root, file)
                    
                    # Create BIDS structure and move the file
                    create_bids_structure(subject_id, session_alt, scan_number, cohort, diagnosis, gender, input_file, path_to_dataset, bids_dir, path_to_clinical)
                    
                    # Write the extracted information to CSV
                    bids_filename = f"sub-MIRIAD{subject_id}_ses-{session}_T1w.nii.gz"
                    output_file = os.path.join(f"sub-MIRIAD{subject_id}", f"ses-{session}", 'anat', bids_filename)
                    csvwriter.writerow([cohort, subject_id, diagnosis, gender, session, input_file, output_file])

                     # Track the minimum age for the participant for baseline
                    if subject_id not in participants_data or participants_data[subject_id]['age'] > age:
                        participants_data[subject_id] = {'participant_id': f"sub-MIRIAD{subject_id}", 
                                                        'sex': gender_clinical, 
                                                        'diagnosis': group, 
                                                        'age': age}

                    # Prepare sessions data
                    sessions_data.append([f"sub-MIRIAD{subject_id}", f"ses-{session}", age])

                    
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
