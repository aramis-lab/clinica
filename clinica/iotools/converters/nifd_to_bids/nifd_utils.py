# coding: utf8

from clinica.iotools.converters.nifd_to_bids.utils.conv_image_folders import (
    dict_conversion,
    get_all_med_name,
)


def break_path(path):
    return path.split("/")


def filter(l_paths, source_dir, descriptors):
    """From a list of paths, remove paths that do not lead to an image described by any of the descriptor in descriptors.

    Args:
        l_paths: list of paths
        source_dir: path to the NIFD dataset
        descriptors: list of descriptor instances

    Returns:
        list of paths
    """
    medical_images = get_all_med_name(source_dir)

    equivalences = dict_conversion(medical_images, descriptors)
    sol = []
    for path in l_paths:
        for elt in break_path(path):
            if elt in equivalences:
                sol.append(path)
                break
    return sol


def csv_to_path(csv_line, source_dir):
    import os

    def change_format_date(date):
        sol = ""
        new_date = date.split("/")
        if len(new_date[0]) != 2:
            new_date[0] = "0" + new_date[0]
        if len(new_date[1]) != 2:
            new_date[1] = "0" + new_date[1]
        sol = new_date[2] + "-" + new_date[0] + "-" + new_date[1]
        return sol

    split = csv_line.split("\t")
    if split[10] != "Original" and split[9] == "MT1; GradWarp; N3m <- Sag 3D MP-RAGE":
        split[9] = "MT1__GradWarp__N3m"
    split[9] = split[9].replace(" ", "_")
    split[9] = split[9].replace(":", "_")
    split[9] = split[9].replace("(", "_")
    split[9] = split[9].replace(")", "_")
    split[9] = split[9].replace("/", "_")
    split[9] = split[9].replace(
        "MT1;_GradWarp;_N3m_<-_Sag_3D_MP-RAGE", "MT1__GradWarp__N3m"
    )
    folder = os.path.join(source_dir, split[0], split[9])

    try:
        subfold = [f.path.split("/")[-1] for f in os.scandir(folder) if f.is_dir()]
    except Exception:
        return None
    date = change_format_date(split[5])

    if (
        split[0] == "1_S_0005"
        and date == "2011-06-13"
        and split[9] in ("NIFD_DTI_b=1000_2.7mm3", "dti_b=2000_64dir", "dti_b=0_scans")
    ):
        date = "2011-06-15"
    if (
        split[0] == "1_S_0084"
        and date == "2012-02-13"
        and split[9]
        in (
            "NIFD_DTI_b1000_2.7mm3",
            "dti_WIP_b2000_64dir",
            "DTI_64_2.2iso_full_ky_fov220",
            "DTI_b=0_2.2iso_full_ky_-_10_acqs",
            "dti_WIP_b0_scans_aah",
        )
    ):
        date = "2012-02-15"
    if (
        split[0] == "1_S_0316"
        and date == "2013-10-30"
        and split[9]
        in (
            "NIFD_DTI_b1000_2.7mm3_511E",
            "ep2d-advdiff-511E_b2000_64dir",
            "ep2d-advdiff-511E_b0_scan",
        )
    ):
        date = "2013-10-31"

    val = [fold for fold in subfold if fold.startswith(date)]

    if val == []:
        return None
    path_date = os.path.join(folder, val[0])
    full_path = os.path.join(
        path_date, [f.path for f in os.scandir(path_date) if f.is_dir()][0]
    )
    return full_path


def get_patients_source_files(source_dir, path_ida):
    """Return a dictionary containing the paths of all Dicom files for a patient.

    Args:
        source_dir: path to the NIFD dataset
        path_ida : path to the ida.tsv file

    Returns:
        patients_source_files[subject_ID] = [paths_to_all_medical_images_of_subject]
    """
    sol = dict()
    fich = open(path_ida, "r")
    fich.readline()
    for line in fich.readlines():
        patient_id = line.split("\t")[0]
        if patient_id not in sol:
            sol[patient_id] = []
        path = csv_to_path(line, source_dir)
        if path is not None:
            sol[patient_id].append(path)
    return sol


def filter_patients_source_files(patients_files, source_dir, descriptors):
    """
    Iterates over a dictionary containing all source files for a patient,
    returns the same dictionary with only path to files that could be converted.

    Args:
        patients_files: patients_source_files[subject_ID] = [paths_to_all_medical_images_of_subject]
        source_dir: path to the NIFD dataset
        descriptors: list of descriptor instances

    Returns:
        patients_source_files[subject_ID] = [paths_to_medical_images_described_by_at_least_one_descriptor]
    """
    for patient in patients_files:
        patients_files[patient] = filter(
            patients_files[patient], source_dir, descriptors
        )

    return patients_files


def create_folder(path):
    """If the provided path does not exist, the required folders are created."""
    import os

    try:
        os.makedirs(path, exist_ok=True)
    except Exception:
        pass


def extract_date(path):
    """Extract the name of the date associated with the given path.

    Args:
        path: Path to a DICOM file

    Returns:
        Date
    """
    path = break_path(path)
    for i in path:
        if i.startswith("20"):
            sol = i.split("_")[0]
            return sol


def extract_name_med_img(path, equivalences):
    """Extract the name of the medical image associated with the given path.

    Args:
        path: Path to a Dicom file
        equivalences: dictionary, equivalences['medical_image_name'] = (Descriptor_instance, modalityLabel)

    Returns:
        Medical name
    """
    path = break_path(path)
    for i in path:
        if i in equivalences:
            return i


def print_patients_source_files(dict):
    """Pretty printer for a nested dictionary such as patients_source_files.

    Args:
        dict: Data structure of the form dict['subject_ID'] = [paths/to/dcm]
    """
    from clinica.utils.stream import cprint

    for key in dict:
        cprint(key)
        for elt in dict[key]:
            cprint("\t" + str(elt))


def print_orderedBIDS(dict):
    """Pretty printer for a nested dictionary such as ordered_bids, final_bids.

    Args:
        dict: Data structure of the form dict['session_id']['dataType']['Priority']['Final_name'] = [paths/to/dcm]
    """
    from clinica.utils.stream import cprint

    for ses in dict:
        cprint(ses)
        for dType in dict[ses]:
            cprint("\t" + str(dType))
            for priority in dict[ses][dType]:
                cprint("\t\t" + str(priority))
                for name in dict[ses][dType][priority]:
                    cprint("\t\t\t" + str(name))
                    for path in dict[ses][dType][priority][name]:
                        cprint("\t\t\t\t" + str(path))


def collect_conversion_tuples(final_bids, dest_dir, patient):
    """Return a list of tuples that will be the input of the convert function.

    If final_bids[session][datatype][priority][name] contains a path, it is added to
    tuple[0], tuple[1] contains the path where the Nifti will be located (follows the BIDS format).

    Args:
        final_bids: A data structure of the form : final_bids[session][datatype][priority][name]
        dest_dir: Home of the BIDS directory
        patient: Instantiation of the Patient class

    Returns:
        A list of tuples [(path/to/dicom, path/to/nifti), ...]
    """
    import os

    sol = []
    for ses in final_bids:
        for dType in final_bids[ses]:
            for priority in final_bids[ses][dType]:
                for name in final_bids[ses][dType][priority]:
                    if final_bids[ses][dType][priority][name] != []:
                        path_source = final_bids[ses][dType][priority][name][0]
                        path_dest = os.path.join(
                            dest_dir, patient.get_name(), ses, dType, name
                        )
                        sol.append((path_source, path_dest))
    return sol


def supress_stdout(func):
    """Wrapper, makes a function non-verbose.

    Args:
        func: function to be silenced
    """
    import contextlib
    import os

    def wrapper(*a, **ka):
        with open(os.devnull, "w") as devnull:
            with contextlib.redirect_stdout(devnull):
                func(*a, **ka)

    return wrapper


def convert_dcm_to_nii(single_tuple):
    """[Summary].

    Args:
        single_tuple: tuple where tuple[0] is the path to the data,
            and tuple[1] the path to the coverted data
    """
    import os
    import subprocess

    from clinica.utils.stream import cprint

    filename = os.path.basename(single_tuple[1])
    path_dest = os.path.dirname(single_tuple[1])
    create_folder(path_dest)
    command = f"dcm2niix -b y -z y -o {path_dest} -f {filename} {single_tuple[0]}"
    cprint(
        msg=f'Converting {os.path.basename(filename).replace("_", " ")}', lvl="debug"
    )
    output_dcm2niix = subprocess.run(command, shell=True, capture_output=True)
    if output_dcm2niix.returncode != 0:
        cprint(
            msg=f'WARNING: errors may have occured.\nOutput:{output_dcm2niix.stdout.decode("utf-8")}\nError:{output_dcm2niix.stderr.decode("utf-8")}',
            lvl="warning",
        )


def convert(list_tuples):
    """Convert a list of tuples = [(path/to/dicom, path/to/nifti), ...].

    Args:
        list_tuples: list of tuples, tuple[0] contains a path to a Dicom file, tuple[1] contains the path where the Nifti file needs to be created.
    """
    from multiprocessing import Pool, cpu_count

    p = Pool(max(cpu_count() - 1, 1))
    p.map(convert_dcm_to_nii, list_tuples)
