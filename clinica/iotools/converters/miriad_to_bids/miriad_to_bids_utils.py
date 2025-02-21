import os
import shutil
import nibabel as nib

# Helper function to create BIDS folders and move files
def create_bids_structure(subject_id, session, run_label, cohort, diagnosis, gender, input_file, path_to_dataset, output_dir, path_to_clinical):
    """Create BIDS folder structure and move files into it."""
    sub_id = f"sub-MIRIAD{subject_id}"
    ses_id = f"ses-{session}"
    run_id = f"run-{run_label}"  # Run number (e.g., run-01)

    # Create output directory for this subject/session
    anat_dir = os.path.join(output_dir, sub_id, ses_id, 'anat')
    os.makedirs(anat_dir, exist_ok=True)
    
    # Convert the input file to .nii.gz if necessary
    input_file_gz = convert_to_nii_gz(input_file)
    
    # Destination filename in BIDS format with run number
    bids_filename = f"{sub_id}_{ses_id}_{run_id}_T1w.nii.gz"
    
    # Copy and rename the file to BIDS format
    shutil.copy(input_file_gz, os.path.join(anat_dir, bids_filename))


# Function to extract subject, session, and run info from filenames
def parse_filename(filename):
    parts = filename.split('_')
    cohort_name = parts[0]  # "miriad"
    subject_id = parts[1]   # e.g., "215"
    diagnosis = parts[2]    # e.g., "AD" or "HC"
    gender = parts[3]       # "M" or "F"
    session_id = parts[4]   # e.g., "01"
    modality = parts[5]     # e.g., "MR"
    run_id = parts[6]       # e.g., "1" (for run-01, run-02)

    return subject_id, session_id, run_id

def convert_to_nii_gz(input_file):
    """Convert a .nii file to .nii.gz format without deleting the original .nii file."""
    if input_file.endswith(".nii.gz"):
        return input_file
    img = nib.load(input_file)
    output_file = input_file.replace(".nii", ".nii.gz")
    nib.save(img, output_file)
    return output_file