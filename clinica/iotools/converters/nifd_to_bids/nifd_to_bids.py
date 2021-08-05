# coding: utf8

"""Convert the NIFD dataset into BIDS."""


def convert_images(path_to_dataset, bids_dir, path_to_clinical):
    """Convert the entire dataset in BIDS.

    Scans available files in the path_to_dataset,
    identifies the patients that have images described by the JSON file,
    converts the image with the highest quality for each category.
    """
    import os

    from clinica.utils.stream import cprint

    from .nifd_utils import (
        collect_conversion_tuples,
        convert,
        filter_patients_source_files,
        get_patients_source_files,
    )
    from .preprocessing.parse_ida import process_ida
    from .preprocessing.update_clinical_info import update_info_clinical
    from .utils.conv_image_folders import (
        dict_conversion,
        get_all_med_name,
        get_descriptors,
    )
    from .utils.manage_conflicts import Manage_conflicts
    from .utils.patient import Patient

    path_converter = os.path.join(os.path.dirname(os.path.realpath(__file__)))
    path_conflicts = os.path.join(
        path_converter, "config_files", "unique_conflicts.txt"
    )
    path_to_clinical_file = os.path.join(
        path_to_clinical, "NIFD_Clinical_Data_2017_final_updated.xlsx"
    )

    path_to_ida = os.path.join(path_to_clinical, "ida.tsv")

    path_to_clinical_info = os.path.join(path_to_clinical, "clinical_info.tsv")

    name_ida = None
    for filename in os.listdir(path_to_clinical):
        if filename.lower().startswith("idasearch") and filename.lower().endswith(
            ".csv"
        ):
            name_ida = filename
    if not name_ida:
        name_ida = "idaSearch_all.csv"

    path_idaSearch = os.path.join(path_to_clinical, name_ida)

    def get_data_dictionary(path_to_clinical_data_folder):
        """Temporary function to get DataDictionary_NIFD_2017.10.18.xlsx file.

        See https://github.com/aramis-lab/clinica/issues/122 for details.

        Args:
            path_to_clinical_data_folder: Path to clinical data folder.

        Returns:
            Path to 'DataDictionary_NIFD_2017.10.18.xlsx' file.
        """
        import os

        from clinica.utils.inputs import RemoteFileStructure, get_file_from_server

        local_nifd_dictionary = os.path.join(
            path_to_clinical, "DataDictionary_NIFD_2017.10.18.xlsx"
        )
        if os.path.exists(local_nifd_dictionary):
            path_to_nifd_dictionary = local_nifd_dictionary
        else:
            NIFD_DICTIONNARY = RemoteFileStructure(
                filename="DataDictionary_NIFD_2017.10.18.xlsx",
                url="https://aramislab.paris.inria.fr/files/data/databases/converters/",
                checksum="e75b23a9f4dad601463f48031cfc00e1180e4877d0bebbdfd340fdbcbacab5cb",
            )
            path_to_nifd_dictionary = get_file_from_server(NIFD_DICTIONNARY)

        return path_to_nifd_dictionary

    path_DataDictionary_NIFD_2017 = get_data_dictionary(path_to_clinical)

    # Pre-processing step, to be executed the first time the converter is used.
    if not os.path.isfile(path_to_ida):
        if os.path.isfile(path_idaSearch):
            cprint(
                "ida.tsv was not found in the clinical data directory, "
                "ida.tsv will be created from idaSearch_<date_of_download>_all.csv"
            )

        else:
            cprint(
                "\nida.tsv does not exist and idaSearch_<date_of_download>_all.csv "
                "was not found in the clinical data directory,"
                " to create it please enter path/to/idaSearch_all.csv :"
            )
            path_idaSearch = input()
            path_idaSearch = path_idaSearch.strip(" ")
        cprint("Creating ida.tsv ...")
        process_ida(path_idaSearch, path_to_clinical)
        assert os.path.isfile(path_to_ida), "Failed to create ida.tsv"
        cprint("ida.tsv successfully created !")

    if not os.path.isfile(path_to_clinical_info):
        if os.path.isfile(path_idaSearch):
            cprint(
                "clinical_info.tsv was not found in the clinical data directory, "
                "clinical_info.tsv will be created from DataDictionary_NIFD_2017.10.18.xlsx"
            )

        else:
            cprint(
                "\nclinical_info.tsv does not exist and DataDictionary_NIFD_2017.10.18.xlsx was not found in the "
                "clinical data directory"
                ", to create it please enter path/to/DataDictionary_NIFD_2017.10.18.xlsx :"
            )
            path_DataDictionary_NIFD_2017 = input()
            path_DataDictionary_NIFD_2017 = path_DataDictionary_NIFD_2017.strip(" ")
        cprint("Creating clinica_info.tsv ...")

        path_clinicals = os.path.join(
            path_converter, "..", "clinical_data_bids_correspondence"
        )

        update_info_clinical(
            path_DataDictionary_NIFD_2017,
            path_clinicals,
            path_to_clinical_file,
            path_to_clinical,
        )
        assert os.path.isfile(
            path_to_clinical_info
        ), "Failed to create clinical_info.tsv"
        cprint("clinical_info.tsv successfully created !")

    cprint("Parsing files to be converted...")
    medical_images = get_all_med_name(path_to_dataset)
    assert len(medical_images) > 0, "The dataset contains no medical image"

    descriptors = get_descriptors(os.path.join(path_converter, "config_files"))
    assert len(descriptors) > 0, "Failed to load the descriptors"

    # equivalences['medical_image_name'] = (Descriptor_instance, modalityLabel)
    # ex : equivalences['T1_mprage_S3_DIS3D'] -> (<Descriptor.Descriptor object at 0x105459ac8>, 'T1w')
    equivalences = dict_conversion(medical_images, descriptors)

    # patients_source_files[subject_ID] = [paths_to_all_medical_images_of_subject]
    # Only contains files that are available in the provided dataset
    patients_source_files = get_patients_source_files(path_to_dataset, path_to_ida)
    patients_source_files = filter_patients_source_files(
        patients_source_files, path_to_dataset, descriptors
    )
    patients_source_files = {
        pat: patients_source_files[pat]
        for pat in patients_source_files
        if patients_source_files[pat] != []
    }

    cf = Manage_conflicts(path_conflicts)

    to_convert = []
    for pat in patients_source_files:
        path_patient = os.path.join(path_to_dataset, pat)

        patient = Patient(pat, path_patient, path_to_ida)
        final_bids = patient.clean_conflicts(
            equivalences, descriptors, patients_source_files[pat], cf, pat
        )

        to_convert.extend(collect_conversion_tuples(final_bids, bids_dir, patient))

    assert to_convert != [], "No Dicom files to convert!"
    cprint("Converting files to Nifti")

    # Converting only images that have not been already converted,
    # the converter does not have to restart from scratch if something fails
    final_convert = []
    for tuple in to_convert:
        if not os.path.isfile(tuple[1] + ".nii.gz"):
            final_convert.append(tuple)

    convert(final_convert)
    return to_convert


def convert_clinical_data(bids_dir, path_to_clinical, to_convert):
    # clinical specifications in BIDS
    import os

    import clinica.iotools.bids_utils as bids
    from clinica.iotools.converters.nifd_to_bids.utils.parse_clinical import (
        Parse_clinical,
    )
    from clinica.utils.stream import cprint

    path_to_ida = os.path.join(path_to_clinical, "ida.tsv")
    assert os.path.isfile(path_to_ida), "Failed to create ida.tsv"

    cprint("Creating clinical data files")

    bids.write_modality_agnostic_files("NIFD", bids_dir)
    pc = Parse_clinical(path_to_clinical)
    pc.make_all(bids_dir)
    pc.make_all_scans(to_convert)
