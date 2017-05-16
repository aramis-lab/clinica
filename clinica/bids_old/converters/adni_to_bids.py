
from clinica.bids.abstract_converter import Converter
from clinica.engine.cmdparser import CmdParser
import logging

__author__ = "Jorge Samper and Sabrina Fontanella"
__copyright__ = "Copyright 2017, The Aramis Lab Team"
__credits__ = [""]
__license__ = ""
__version__ = "1.0.0"
__maintainer__ = "Sabrina Fontanella"
__email__ = "sabrina.fontanella@icm-institute.org"
__status__ = "Development"


class ADNI_TO_BIDS(Converter, CmdParser):

    def define_name(self):
        self._name = 'adni-to-bids'

    def define_options(self):
        self._args.add_argument("dataset_directory",
                               help='Path of the unorganized ADNI directory.')
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')
        self._args.add_argument("--clinical-only", type=bool, default=False,
                                dest='clinical_only',
                                help='(Optional) Given an already existing BIDS output folder, convert only the clinical data.')

    def convert_clinical_data(self, input_path, out_path):
        """

        Convert clinical data of ADNI dataset.

        :param src:
        :param out_path:
        :return:

        {'sub-ADNI009S1354': { 'bl': {u'diagnosis': 'Dementia', 'session_id': 'ses-bl'},
                              'm06': {u'diagnosis': 'Dementia', 'session_id': 'ses-m06'}},...
         ....
         }
        """
        from os import path
        import os
        import logging
        import clinica.bids.bids_utils as bids
        import pandas as pd
        from glob import glob
        from os.path import normpath

        logging.basicConfig(filename=path.join(out_path, 'conversion_clinical.log'),
                            format='%(asctime)s %(levelname)s:%(message)s', level=logging.DEBUG,
                            datefmt='%m/%d/%Y %I:%M')

        clinic_specs_path = path.join(os.path.dirname(os.path.dirname(__file__)), 'data',
                                      'clinical_specifications_adni.xlsx')


        try:
            os.path.exists(out_path)
        except IOError:
            print 'BIDS folder not found.'
            raise

        bids_ids = bids.get_bids_subjs_list(out_path)
        bids_subjs_paths = bids.get_bids_subjs_paths(out_path)

        # -- Creation of participant.tsv --
        participants_df = bids.create_participants_df(input_path, out_path, 'ADNI', 'clinical_specifications_adni.xlsx', bids_ids)
        participants_df.to_csv(path.join(out_path, 'participant.tsv'), sep='\t', index=False)

        # -- Creation of sessions.tsv --
        logging.info("--Creation of sessions files. --")
        print("\nCreation of sessions files...")
        # Load data
        sessions = pd.read_excel(clinic_specs_path, sheetname='sessions.tsv')
        sessions_fields = sessions['ADNI']
        field_location = sessions['ADNI location']
        sessions_fields_bids = sessions['BIDS CLINICA']
        file_to_read = pd.read_csv(path.join(input_path, 'clinicalData', 'ADNIMERGE.csv'))
        fields_dataset = []
        fields_bids = []
        sessions_dict = {}
        subj_to_remove = []

        for i in range(0, len(sessions_fields)):
            if not pd.isnull(sessions_fields[i]):
                fields_bids.append(sessions_fields_bids[i])
                fields_dataset.append(sessions_fields[i])

        sessions_df = pd.DataFrame(columns=fields_bids)

        # For each line of ADNIMERGE check if the subjects is available in the BIDS dataset and extract the information
        for r in range(0, len(file_to_read.values)):
            row = file_to_read.iloc[r]
            id_ref = 'PTID'
            subj_id = row[id_ref.decode('utf-8')]
            # Removes all the - from
            subj_id_alpha = bids.remove_space_and_symbols(subj_id)
            subj_bids = [s for s in bids_ids if subj_id_alpha in s]
            if len(subj_bids) == 0:
                pass
            else:
                subj_bids = subj_bids[0]
                for i in range(0, len(sessions_fields)):
                    # If the i-th field is available
                    if not pd.isnull(sessions_fields[i]):
                        sessions_df[sessions_fields_bids[i]] = row[sessions_fields[i]]
                        visit_id = row['VISCODE']
                        # If the dictionary already contain the subject add or update information regarding a specific session,
                        # otherwise create the enty
                        if sessions_dict.has_key(subj_bids):
                            sess_available = sessions_dict[subj_bids].keys()
                            if visit_id in sess_available:
                                sessions_dict[subj_bids][visit_id].update({sessions_fields_bids[i]: row[sessions_fields[i]]})
                            else:
                                sessions_dict[subj_bids].update({visit_id: {'session_id': 'ses-'+visit_id,
                                                                            sessions_fields_bids[i]: row[sessions_fields[i]]}})
                        else:
                            sessions_dict.update({subj_bids: {visit_id: {'session_id': 'ses-'+visit_id,
                                                                         sessions_fields_bids[i]: row[sessions_fields[i]]}}})

        # Write the information contained inside the dictionary to the proper session file
        for sp in bids_subjs_paths:
            bids_id = sp.split(os.sep)[-1]
            sessions_df = pd.DataFrame(columns=fields_bids)
            if sessions_dict.has_key(bids_id):
                sess_aval = sessions_dict[bids_id].keys()
                for ses in sess_aval:
                    sessions_df = sessions_df.append(pd.DataFrame(sessions_dict[bids_id][ses], index=['i', ]))

                sessions_df.to_csv(path.join(sp, bids_id + '_sessions.tsv'), sep='\t', index=False)

        # -- Creation of scans files --
        print 'Creation of scans files...'
        scans_dict = {}

        for bids_id in bids_ids:
            scans_dict.update({bids_id: {'T1/DWI/fMRI': {}, 'FDG': {}}})

        scans_specs = pd.read_excel(clinic_specs_path, sheetname='scans.tsv')
        scans_fields_db = scans_specs['ADNI']
        scans_fields_bids = scans_specs['BIDS CLINICA']
        scans_fields_mod = scans_specs['Modalities related']
        fields_bids = ['filename']

        for i in range(0, len(scans_fields_db)):
            if not pd.isnull(scans_fields_db[i]):
                fields_bids.append(scans_fields_bids[i])

        scans_df = pd.DataFrame(columns=(fields_bids))

        for bids_subj_path in bids_subjs_paths:
            # Create the file
            bids_id = os.path.basename(normpath(bids_subj_path))

            sessions_paths = glob(path.join(bids_subj_path, 'ses-*'))
            for session_path in sessions_paths:
                session_name = session_path.split(os.sep)[-1]
                tsv_name = bids_id + '_' + session_name + "_scans.tsv"

                # If the file already exists, remove it
                if os.path.exists(path.join(session_path, tsv_name)):
                    os.remove(path.join(session_path, tsv_name))

                scans_tsv = open(path.join(session_path, tsv_name), 'a')
                scans_df.to_csv(scans_tsv, sep='\t', index=False)

                # Extract modalities available for each subject
                mod_available = glob(path.join(session_path, '*'))
                for mod in mod_available:
                    mod_name = os.path.basename(mod)
                    files = glob(path.join(mod, '*'))
                    for file in files:
                        file_name = os.path.basename(file)
                        if mod == "anat" or mod == "dwi" or mod == "func":
                            type_mod = 'T1/DWI/fMRI'
                        else:
                            type_mod = 'FDG'

                        scans_df['filename'] = pd.Series(path.join(mod_name, file_name))
                        scans_df.to_csv(scans_tsv, header=False, sep='\t', index=False)

                scans_df = pd.DataFrame(columns=(fields_bids))

        print '-- Scans files created for each subject. --'

    def convert_images(self, source_dir, clinical_dir, dest_dir, subjs_list_path='', mod_to_add='', mod_to_update='', compute_paths = True):
        """
        The function first computes the paths of the right image to be converted and
        :param source_dir:
        :param dest_dir:
        :return:
        """

        import clinica.bids.bids_utils as bids
        import os
        from os import path
        import pandas as pd
        import logging
        from datetime import datetime
        import adni_utils
        import adni_modalities.adni_dwi as adni_dwi
        import adni_modalities.adni_fmri as adni_fmri

        print "*******************************"
        print "ADNI to BIDS converter"

        if mod_to_update != '':
            print 'Updating: ', mod_to_update



        bids_ids = []
        alpha_ids = []
        adni_merge_path = path.join(clinical_dir, 'ADNIMERGE.csv')
        adni_merge = pd.read_csv(adni_merge_path)


        # Options to use
        convert_func = (mod_to_add == '' or mod_to_add == 'func') and (mod_to_update == '' or mod_to_update == 'func')
        convert_dwi = (mod_to_add == '' or mod_to_add == 'dwi') and (mod_to_update == '' or mod_to_update == 'dwi')
        convert_anat = (mod_to_add == '' or mod_to_add == 'anat') and (mod_to_update == '' or mod_to_update == 'anat')
        convert_pet = (mod_to_add == '' or mod_to_add == 'pet') and (mod_to_update == '' or mod_to_update == 'pet')

        # Load a file with subjects list or compute all the subjects
        if subjs_list_path != '':
            subjs_list = [line.rstrip('\n') for line in open(subjs_list_path)]
        else:
            print 'Using all the subjects contained into ADNI merge...'
            subjs_list = adni_merge['PTID'].unique()

        print 'Subjects found:', len(subjs_list)
        print "*******************************"

        if (mod_to_add == '' and mod_to_update == '') or (mod_to_add != '' and os.path.exists(dest_dir) == False):
            print 'Creating the output folder'
            os.mkdir(dest_dir)
            os.mkdir(path.join(dest_dir, 'conversion_info'))
            os.mkdir(path.join(dest_dir, 'conversion_info', 'log'))

        logging.basicConfig(filename=path.join(dest_dir, path.join(dest_dir, 'conversion_info', 'log', datetime.now().strftime('log_%d_%m_%Y_%H_%M.log'))), format='%(asctime)s %(levelname)s:%(message)s',
                            datefmt='%m/%d/%Y %I:%M', level=logging.DEBUG)

        # Compute anat paths
        if convert_anat:
            t1_paths = adni_utils.compute_t1_paths(source_dir, clinical_dir, dest_dir, subjs_list)

        # Compute pet paths
        if convert_pet:
            print 'Calculating paths for PET FDG...'
            logging.info('Calculating paths for PET FDG...')
            pet_fdg_paths = adni_utils.compute_fdg_pet_paths(source_dir, clinical_dir, dest_dir, subjs_list)
            print 'Done!'
            logging.info('Done')
            pet_av45_paths = adni_utils.compute_av45_pet_paths(source_dir, clinical_dir, dest_dir, subjs_list)

        # Compute func paths
        if convert_func:
            print 'Calculating paths for FMRI...'
            logging.info('Calculating paths for FMRI...')
            #fmri_paths = adni_utils.compute_fmri_path(source_dir, clinical_dir, dest_dir, subjs_list)
            fmri_paths = pd.read_csv(path.join(dest_dir,'conversion_info','fmri_paths.tsv'), sep='\t')
            print 'Done!'
            logging.info('Done!')

        # Compute dwi paths
        if convert_dwi:
            if compute_paths:
                dwi_paths = adni_dwi.compute_dti_paths(source_dir, clinical_dir, dest_dir, subjs_list)
            else:
                print 'Loading the old dwi_paths file...'
                dwi_paths = pd.read_csv(path.join(dest_dir, 'conversion_info', 'dwi_paths.tsv'), sep='\t')

        # Create subjects folders
        for subj in subjs_list:
            alpha_id = bids.remove_space_and_symbols(subj)
            bids_id = 'sub-ADNI' + alpha_id
            alpha_ids.append(alpha_id)
            bids_ids.append(bids_id)
            if (mod_to_add == '' and mod_to_update == '') or (
                    mod_to_add != '' and not os.path.exists(path.join(dest_dir, bids_id))):
                os.mkdir(path.join(dest_dir, bids_id))

        if convert_func:
            if mod_to_update == 'func':
                adni_fmri.convert_fmri(dest_dir, subjs_list, fmri_paths, mod_to_add=False, mod_to_update=True)

        if convert_dwi:
            if mod_to_update == 'dwi':
                adni_dwi.convert_dwi(dest_dir, dwi_paths, mod_to_add=False, mod_to_update=True)
            elif mod_to_add == 'dwi':
                adni_dwi.convert_dwi(dest_dir, dwi_paths, mod_to_add=True, mod_to_update=False)


        print '\n'

    def visits_to_timepoints_dti(self, subject, ida_meta_subj, adnimerge_subj):
        from datetime import datetime

        visits = dict()
        unique_visits = list(ida_meta_subj.Visit.unique())
        pending_timepoints = []

        # We try to obtain the corresponding image Visit for a given VISCODE
        for adni_row in adnimerge_subj.iterrows():
            visit = adni_row[1]
            if visit.ORIGPROT == 'ADNI2':
                if visit.VISCODE == 'bl':
                    preferred_visit_name = 'ADNI2 Screening MRI-New Pt'
                elif visit.VISCODE == 'm03':
                    preferred_visit_name = 'ADNI2 Month 3 MRI-New Pt'
                elif visit.VISCODE == 'm06':
                    preferred_visit_name = 'ADNI2 Month 6-New Pt'
                else:
                    year = str(int(visit.VISCODE[1:]) / 12)
                    preferred_visit_name = 'ADNI2 Year ' + year + ' Visit'
            else:
                if visit.VISCODE == 'bl':
                    if visit.ORIGPROT == 'ADNI1':
                        preferred_visit_name = 'ADNI Screening'
                    else:  # ADNIGO
                        preferred_visit_name = 'ADNIGO Screening MRI'
                elif visit.VISCODE == 'm03':  # Only for ADNIGO Month 3
                    preferred_visit_name = 'ADNIGO Month 3 MRI'
                else:
                    month = int(visit.VISCODE[1:])
                    if month < 54:
                        preferred_visit_name = 'ADNI1/GO Month ' + str(month)
                    else:
                        preferred_visit_name = 'ADNIGO Month ' + str(month)

            if preferred_visit_name in unique_visits:
                key_preferred_visit = (visit.VISCODE, visit.COLPROT, visit.ORIGPROT)
                if key_preferred_visit not in visits.keys():
                    visits[key_preferred_visit] = preferred_visit_name
                elif visits[key_preferred_visit] != preferred_visit_name:
                    print 'Multiple visits for one timepoint!'
                    print subject
                    print key_preferred_visit
                    print visits[key_preferred_visit]
                    print visit
                unique_visits.remove(preferred_visit_name)
                continue

            pending_timepoints.append(visit)

        # Then for images.Visit non matching the expected labels we find the closest date in visits list
        for visit in unique_visits:
            image = (ida_meta_subj[ida_meta_subj.Visit == visit]).iloc[0]
            min_db = 100000
            min_db2 = 0
            min_visit = None
            min_visit2 = None

            for timepoint in pending_timepoints:
                db = self.days_between(image['Scan Date'], timepoint.EXAMDATE)
                if db < min_db:
                    min_db2 = min_db
                    min_visit2 = min_visit

                    min_db = db
                    min_visit = timepoint

            if min_visit is None:
                print 'No corresponding timepoint in ADNIMERGE for subject ' + subject + ' in visit ' + image.Visit
                print image
                continue

            if min_visit2 is not None and min_db > 90:
                print 'More than 60 days for corresponding timepoint in ADNIMERGE for subject ' + subject + ' in visit ' + image.Visit + ' on ' + image.ScanDate
                print 'Timepoint 1: ' + min_visit.VISCODE + ' - ' + min_visit.ORIGPROT + ' on ' + min_visit.EXAMDATE + ' (Distance: ' + str(
                    min_db) + ' days)'
                print 'Timepoint 2: ' + min_visit2.VISCODE + ' - ' + min_visit2.ORIGPROT + ' on ' + min_visit2.EXAMDATE + ' (Distance: ' + str(
                    min_db2) + ' days)'

                # If image is too close to the date between two visits we prefer the earlier visit
                if (datetime.strptime(min_visit.EXAMDATE, "%Y-%m-%d")
                        > datetime.strptime(image.ScanDate, "%Y-%m-%d")
                        > datetime.strptime(min_visit2.EXAMDATE, "%Y-%m-%d")):
                    dif = self.days_between(min_visit.EXAMDATE, min_visit2.EXAMDATE)
                    if abs((dif / 2.0) - min_db) < 30:
                        min_visit = min_visit2

                print 'We prefer ' + min_visit.VISCODE

            key_min_visit = (min_visit.VISCODE, min_visit.COLPROT, min_visit.ORIGPROT)
            if key_min_visit not in visits.keys():
                visits[key_min_visit] = image.Visit
            elif visits[key_min_visit] != image.Visit:
                print 'Multiple visits for one timepoint!'
                print subject
                print key_min_visit
                print visits[key_min_visit]
                print image.Visit

        return visits

    def run_pipeline(self, args):
        if args.modality is True:
            self.convert_clinical_data(args.dataset_directory, args.bids_directory)

    # def dti_to_bids(self, dti_paths, bids_dir):
    #     from clinica.bids.bids_utils import remove_space_and_symbols, dcm_to_nii
    #     import pandas as pd
    #     from os import path, makedirs
    #     from numpy import nan
    #
    #     dti_images = pd.io.parsers.read_csv(dti_paths, sep='\t')
    #
    #     for row in dti_images.iterrows():
    #         image = row[1]
    #         if image.Path is nan:
    #             continue
    #
    #         subject = 'sub-' + remove_space_and_symbols(image.Subject_ID)
    #         session = 'ses-' + self.viscode_to_session(image.VISCODE)
    #
    #         output_path = path.join(bids_dir, subject, session, 'dwi')
    #         #TODO Define the standard notation for acquisition name: axial and axialEnhanced? Or enhancedAxial?
    #         bids_name = subject + '_' + session + '_acq-' + ('axialEnhanced' if image.Enhanced else 'axial') + '_dwi'
    #
    #         try:
    #             makedirs(output_path)
    #         except OSError:
    #             if not path.isdir(output_path):
    #                 raise
    #         dcm_to_nii(image.Path, output_path, bids_name)
