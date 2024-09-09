"""Convert MIRIAD dataset to BIDS."""

from pathlib import Path
from typing import Optional

import os
import shutil
import csv
from clinica.utils.filemanip import UserProvidedPath

# Paths
input_dir = 'your_input_directory'  # Where the original data is located
output_dir = 'your_output_directory'  # Where the BIDS data will be written
csv_file = 'metadata.csv'  # Metadata CSV file to store extracted information

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
    from clinica.iotools.converters.miriad_to_bids.miriad_to_bids_utils import create_bids_structure
    # Prepare CSV
    with open(csv_file, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['cohort', 'subject_id', 'diagnosis', 'gender', 'session', 'input_file', 'output_file'])
        
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
                    session = parts[4]    # Session number
                    
                    # Full path of input file
                    input_file = os.path.join(root, file)
                    
                    # Create BIDS structure and move the file
                    create_bids_structure(subject_id, session, cohort, diagnosis, gender, input_file, path_to_dataset, bids_dir, path_to_clinical)
                    
                    # Write the extracted information to CSV
                    bids_filename = f"sub-{subject_id}_ses-{session}_T1w.nii.gz"
                    output_file = os.path.join(f"sub-{subject_id}", f"ses-{session}", 'anat', bids_filename)
                    csvwriter.writerow([cohort, subject_id, diagnosis, gender, session, input_file, output_file])

    print("Conversion to BIDS format and metadata extraction completed.")
