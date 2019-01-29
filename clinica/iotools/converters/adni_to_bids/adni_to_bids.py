# coding: utf-8

from clinica.iotools.abstract_converter import Converter

__author__ = "Jorge Samper-Gonzalez and Sabrina Fontanella"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = [""]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper-Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


class AdniToBids(Converter):

    def get_modalities_supported(self):
        """
        Return a list of modalities supported

        Returns: a list containing the modalities supported by the converter (T1, PET_FDG, PET_AV45)

        """
        return ['T1', 'PET_FDG', 'PET_AV45', 'DWI', 'fMRI', 'FLAIR']

    def check_adni_dependencies(self):
        """
        Check the dependencies of ADNI converter
        """
        from clinica.utils.check_dependency import is_binary_present

        list_binaries = ['dcm2nii', 'dcm2niix']

        for binary in list_binaries:
            if not is_binary_present(binary):
                raise RuntimeError(
                    '%s is not present '
                    'in your PATH environment.' % binary)

    def convert_clinical_data(self,  clinical_data_dir, out_path):
        """
        Convert the clinical data of ADNI specified into the file clinical_specifications_adni.xlsx

        Args:
            clinical_data_dir:  path to the clinical data directory
            out_path: path to the BIDS directory

        """
        from os import path
        import os
        import pandas as pd
        import clinica.iotools.bids_utils as bids
        from clinica.utils.stream import cprint
        import clinica.iotools.converters.adni_to_bids.adni_utils as adni_utils

        clinic_specs_path = path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'data',
                                      'clinical_specifications_adni.xlsx')
        try:
            os.path.exists(out_path)
        except IOError:
            print('BIDS folder not found.')
            raise

        bids_ids = bids.get_bids_subjs_list(out_path)
        bids_subjs_paths = bids.get_bids_subjs_paths(out_path)

        if not bids_ids:
            adni_merge_path = path.join(clinical_data_dir, 'ADNIMERGE.csv')
            adni_merge = pd.io.parsers.read_csv(adni_merge_path, sep=',')
            bids_ids = ['sub-ADNI' + subj.replace('_', '') for subj in list(adni_merge.PTID.unique())]
            bids_subjs_paths = [path.join(out_path, subj) for subj in bids_ids]

        # -- Creation of participant.tsv --
        cprint("Creating participants.tsv...")
        participants_df = bids.create_participants_df('ADNI', clinic_specs_path, clinical_data_dir, bids_ids)

        # Replace the original values with the standard defined by the AramisTeam
        participants_df['sex'] = participants_df['sex'].replace('Male', 'M')
        participants_df['sex'] = participants_df['sex'].replace('Female', 'F')

        participants_df.to_csv(path.join(out_path, 'participants.tsv'), sep='\t', index=False, encoding='utf-8')

        # -- Creation of sessions.tsv --
        cprint("Creating sessions files...")
        adni_utils.create_adni_sessions_dict(bids_ids, clinic_specs_path, clinical_data_dir, bids_subjs_paths)

        # -- Creation of scans files --
        cprint('Creating scans files...')
        adni_utils.create_adni_scans_files(clinic_specs_path, bids_subjs_paths, bids_ids)

    def convert_images(self, source_dir, clinical_dir, dest_dir, subjs_list_path='', modalities=['T1', 'PET_FDG', 'PET_AV45', 'DWI']):
        """
        Convert the images of ADNI

        Args:
            source_dir: path to the ADNI directory
            clinical_dir: path to the clinical data directory
            dest_dir: path to the BIDS directory
            subjs_list_path: list of subjects to process
            modalities: modalities to convert (T1, PET_FDG, PET_AV45, DWI, FLAIR)

        """

        import os
        from os import path
        import pandas as pd
        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_t1 as adni_t1
        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_av45_pet as adni_av45
        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fdg_pet as adni_fdg
        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_dwi as adni_dwi
        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fmri as adni_fmri
        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_flair as adni_flair
        from clinica.utils.stream import cprint

        adni_merge_path = path.join(clinical_dir, 'ADNIMERGE.csv')
        adni_merge = pd.read_csv(adni_merge_path)

        # Load a file with subjects list or compute all the subjects
        if subjs_list_path is not None and subjs_list_path != '':
            cprint('Loading a subjects lists provided by the user...')
            subjs_list = [line.rstrip('\n') for line in open(subjs_list_path)]
        else:
            cprint('Using all the subjects contained into the ADNIMERGE.csv file...')
            subjs_list = list(adni_merge['PTID'].unique())

        # Create the output folder if is not already existing
        if not os.path.exists(dest_dir):
            os.makedirs(dest_dir)
            os.makedirs(path.join(dest_dir, 'conversion_info'))

        if 'T1' in modalities:
            adni_t1.convert_adni_t1(source_dir, clinical_dir, dest_dir, subjs_list)

        if 'PET_FDG' in modalities:
            adni_fdg.convert_adni_fdg_pet(source_dir, clinical_dir, dest_dir, subjs_list)

        if 'PET_AV45' in modalities:
            adni_av45.convert_adni_av45_pet(source_dir, clinical_dir, dest_dir, subjs_list)

        if 'DWI' in modalities:
            adni_dwi.convert_adni_dwi(source_dir, clinical_dir, dest_dir, subjs_list)

        if 'fMRI' in modalities:
            adni_fmri.convert_adni_fmri(source_dir, clinical_dir, dest_dir, subjs_list)

        if 'FLAIR' in modalities:
            adni_flair.convert_adni_flair(source_dir, clinical_dir, dest_dir, subjs_list)
