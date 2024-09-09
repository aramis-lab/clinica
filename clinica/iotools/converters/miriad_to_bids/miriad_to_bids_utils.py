import os
import shutil

# Helper function to create BIDS folders and move files
def create_bids_structure(subject_id, session, cohort, diagnosis, gender, input_file, path_to_dataset, output_dir, path_to_clinical
):

    """_summary_

    Args:
        session (_type_): _description_
        cohort (_type_): _description_
        diagnosis (_type_): _description_
        gender (_type_): _description_
        input_file (_type_): _description_
        output_dir (_type_): _description_
        path_to_dataset (_type_, optional): _description_. Defaults to None, n_procs: Optional[int] = 1, **kwargs, ):#subject_id.
    """
    sub_id = f"sub-{subject_id}"
    ses_id = f"ses-{session}"
    
    # Create output directory for this subject/session
    anat_dir = os.path.join(output_dir, sub_id, ses_id, 'anat')
    os.makedirs(anat_dir, exist_ok=True)
    
    # Destination filename in BIDS format
    bids_filename = f"{sub_id}_{ses_id}_T1w.nii.gz"
    
    # Copy and rename the file to BIDS format
    shutil.copy(input_file, os.path.join(anat_dir, bids_filename))
