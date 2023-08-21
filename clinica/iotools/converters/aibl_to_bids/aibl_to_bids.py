"""
Convert the AIBL dataset (https://www.aibl.csiro.au/) into BIDS.
"""


def convert_images(
    path_to_dataset: str, path_to_csv: str, bids_dir: str, overwrite: bool = False
) -> None:
    """Conversion of the entire dataset in BIDS."""
    from os.path import exists
    from pathlib import Path

    from clinica.iotools.converters.aibl_to_bids.utils import Modality, paths_to_bids
    from clinica.utils.stream import cprint

    path_to_dataset = Path(path_to_dataset)
    path_to_csv = Path(path_to_csv)
    bids_dir = Path(bids_dir)

    list_of_created_files = [
        paths_to_bids(
            path_to_dataset,
            path_to_csv,
            bids_dir,
            modality,
            overwrite=overwrite,
        )
        for modality in Modality
    ]
    missing_files = []
    for modality_list in list_of_created_files:
        for file in modality_list:
            if not exists(str(file)):
                missing_files.append(file)
    if missing_files:
        msg = "The following file were not converted:\n" + "\n".join(missing_files)
        cprint(msg=msg, lvl="warning")


def convert_clinical_data(bids_dir: str, path_to_csv: str) -> None:
    # clinical specifications in BIDS
    from os.path import join, realpath, split

    import clinica.iotools.bids_utils as bids
    from clinica.iotools.converters.aibl_to_bids.utils import (
        create_participants_tsv_file,
        create_scans_tsv_file,
        create_sessions_tsv_file,
    )
    from clinica.utils.stream import cprint

    clinical_spec_path = join(
        split(realpath(__file__))[0], "../../data/clinical_specifications"
    )
    # if not exists(clinical_spec_path):
    #    raise FileNotFoundError(
    #        f"{clinical_spec_path} file not found ! This is an internal file of Clinica."
    #    )

    cprint("Creating modality agnostic files...")
    readme_data = {
        "link": "http://adni.loni.usc.edu/study-design/collaborative-studies/aibl/",
        "desc": (
            "The Australian Imaging, Biomarker & Lifestyle Flagship Study of Ageing (AIBL) seeks to discover which "
            "biomarkers, cognitive characteristics, and health and lifestyle factors determine the development of AD. "
            "Although AIBL and ADNI have many of the same goals, there are differences between the two projects."
        ),
    }
    bids.write_modality_agnostic_files(
        study_name="AIBL", readme_data=readme_data, bids_dir=bids_dir
    )

    cprint("Creating participants.tsv...")
    create_participants_tsv_file(
        bids_dir, clinical_spec_path, path_to_csv, delete_non_bids_info=True
    )

    cprint("Creating sessions files...")
    create_sessions_tsv_file(bids_dir, path_to_csv, clinical_spec_path)

    cprint("Creating scans files...")
    create_scans_tsv_file(bids_dir, path_to_csv, clinical_spec_path)
