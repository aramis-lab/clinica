# coding: utf8


class Patient(object):
    """Class that handles sessions ordering, image quality hierarchy and final BIDS structure."""

    def __init__(self, name, path_patient, path_ida):
        self.name = name
        self.path = path_patient
        self.path_ida = path_ida

        self.sessions, self.corr_sessions = self.get_sessions()
        self.ses0 = self.get_ses0()

    def get_sessions(self):
        """Return the list of all sessions that the patient attended."""
        import pandas as pd

        def change_format_date(date):
            sol = ""
            new_date = date.split("/")
            if len(new_date[0]) != 2:
                new_date[0] = "0" + new_date[0]
            if len(new_date[1]) != 2:
                new_date[1] = "0" + new_date[1]
            sol = new_date[2] + "-" + new_date[0] + "-" + new_date[1]
            return sol

        df = pd.read_csv(self.path_ida, sep="\t")
        patient_clinical = df.loc[df["Subject ID"] == self.name]
        dic = dict()
        sessions = list(set(list(patient_clinical["Study Date"])))

        for ses in sessions:
            visit = patient_clinical.loc[patient_clinical["Study Date"] == ses][
                "Visit"
            ].iloc[0]
            ses_number = visit.split(" ")[1]
            dic[change_format_date(ses)] = int(ses_number)

        return sessions, dic

    def get_ses0(self):
        if self.sessions == []:
            print("Warning: sessions not defined for " + str(self.name))
            return None
        return self.sessions[0]

    def get_sesNumber(self, ses_date):
        return self.corr_sessions[ses_date]

    def get_sesName(self, ses_date):
        number = str(self.get_sesNumber(ses_date))
        while len(number) < 2:
            number = "0" + number

        res = "ses-M" + number

        return res

    def get_name(self):
        return "sub-NIFD" + self.name.replace("_", "")

    def order_sessions(self, equivalences, descriptors, folders):
        """Order sessions of a given patient.

        Args:
            folders : List of all paths to medical images

        Returns:
            same_ses : dictionary containing paths ordered by sessions
                i.e. same_ses = {'session_number' : [all paths to MRIs made during said session]}
        """
        from clinica.iotools.converters.nifd_to_bids.nifd_utils import extract_date

        dates = [extract_date(path) for path in folders]
        dates = list(set(dates))

        same_dates = [
            (self.get_sesName(date), [fold for fold in folders if date in fold])
            for date in dates
        ]

        ses_names = [tupl[0] for tupl in same_dates]
        ses_names = list(set(ses_names))

        same_ses = {}
        for tupl in same_dates:
            if tupl[0] not in same_ses:
                same_ses[tupl[0]] = tupl[1]
            else:
                same_ses[tupl[0]].extend(tupl[1])

        return same_ses

    def order_priorities(self, equivalences, descriptors, folders):
        """Order all folders in the BIDS format following the priorities defined in the JSON file.

        Args:
            equivalences:   Data structure of the form :
                            equivalences['medical_image_name'] = (Descriptor_instance, modalityLabel),
                            links a medical name to its descriptor
            descriptors: List of descriptor objects, this list is build from the config_dcm2bids.json file
            folders: List of directories' paths containing Dicom files for a given patient

        Returns:
          sol: dictionary of the form sol['session_id']['dataType']['Priority']['Final_name'] = [paths/to/dcm]
        """
        from clinica.iotools.converters.nifd_to_bids.nifd_utils import (
            extract_name_med_img,
        )

        sol = {}
        sessions = self.order_sessions(equivalences, descriptors, folders)

        for ses_id in sessions:
            sol[ses_id] = {}
            for path in sessions[ses_id]:

                med_name = extract_name_med_img(path, equivalences)
                dataType = equivalences[med_name][0].dataType
                priority = equivalences[med_name][0].priority
                final_name = f"{self.get_name()}_{ses_id}_{equivalences[med_name][0].get_bids_info()}"

                if dataType not in sol[ses_id]:
                    sol[ses_id][dataType] = {}
                if priority not in sol[ses_id][dataType]:
                    sol[ses_id][dataType][priority] = {}
                if final_name not in sol[ses_id][dataType][priority]:
                    sol[ses_id][dataType][priority][final_name] = [path]
                elif path not in sol[ses_id][dataType][priority][final_name]:
                    sol[ses_id][dataType][priority][final_name].extend([path])

        return sol

    def clean_conflicts(
        self, equivalences, descriptors, folders, conflicts_manager, pat
    ):
        """Order all folders in the BIDS format and removes all conflicts.

        The output contains all paths to be used by the converter.

        Args:
            equivalences:   Data structure of the form :
                            equivalences['medical_image_name'] = (Descriptor_instance, modalityLabel),
                            links a medical name to its descriptor
            descriptors: List of descriptor objects, this list is build from the config_dcm2bids.json file
            folders: List of directories' paths containing Dicom files for a given patient
            conflicts_manager: Instance of Manage_conflicts(), handles quality of converted images (ex : T1_DIS3D > T1)
            pat: Patient ID

        Returns:
            A data structure of the form : ordered_bids[session][datatype][priority][name]
            ordered_bids[session][datatype][priority][name] ontains a list that is either empty or has a single path to a dicom file.
            All paths will then be converted to Nifti following the BIDS format.
        """
        from clinica.iotools.converters.nifd_to_bids.nifd_utils import (
            extract_name_med_img,
        )

        ordered_bids = self.order_priorities(equivalences, descriptors, folders)

        for ses in ordered_bids:
            for dType in ordered_bids[ses]:
                encountered = []
                highest_priority = max(ordered_bids[ses][dType].keys())
                while highest_priority != 0:
                    if highest_priority in ordered_bids[ses][dType]:
                        for name in ordered_bids[ses][dType][highest_priority]:
                            enc = name.split("_")[-1]
                            if enc in encountered:
                                ordered_bids[ses][dType][highest_priority][name] = []
                            else:
                                encountered.append(enc)
                                if (
                                    len(
                                        ordered_bids[ses][dType][highest_priority][name]
                                    )
                                    > 1
                                ):

                                    conflit = [
                                        extract_name_med_img(path, equivalences)
                                        for path in ordered_bids[ses][dType][
                                            highest_priority
                                        ][name]
                                    ]
                                    try:
                                        select = conflicts_manager.make_decision(
                                            conflit
                                        )
                                    except Exception:
                                        print(
                                            f"Warning : {str(conflit)} not in expected conflicts, files will not be converted [Subject ID: {pat}]"
                                        )
                                        select = None
                                    if not select:
                                        ordered_bids[ses][dType][highest_priority][
                                            name
                                        ] = []
                                    else:
                                        val = [
                                            path
                                            for path in ordered_bids[ses][dType][
                                                highest_priority
                                            ][name]
                                            if str(select) in path
                                        ]
                                        ordered_bids[ses][dType][highest_priority][
                                            name
                                        ] = val
                    highest_priority -= 1

        return ordered_bids

    def get_conflicts(self, equivalences, descriptors, folders, conflicts={}):
        # from nifd_utils import extract_name_med_img
        """
        Returns a list of all possible conflicts in dataset according to the json file
        A conflict happens when several images end up in the same directory with the same name after BIDS conversion
        """
        ordered_bids = self.order_priorities(equivalences, descriptors, folders)

        for ses in ordered_bids:
            for dataType in ordered_bids[ses]:
                highest_priority = max(ordered_bids[ses][dataType].keys())
                name_encountered = []
                while highest_priority != 0:
                    if highest_priority in ordered_bids[ses][dataType]:
                        for name in ordered_bids[ses][dataType][highest_priority]:
                            name2 = name
                            name = name.split("_")[-1]
                            if name not in name_encountered:
                                name_encountered.append(name)
                                if (
                                    len(
                                        ordered_bids[ses][dataType][highest_priority][
                                            name2
                                        ]
                                    )
                                    > 1
                                ):
                                    if name not in conflicts:
                                        conflicts[name] = [
                                            ordered_bids[ses][dataType][
                                                highest_priority
                                            ][name2]
                                        ]
                                    else:
                                        conflicts[name].append(
                                            ordered_bids[ses][dataType][
                                                highest_priority
                                            ][name2]
                                        )

                    highest_priority -= 1

        return conflicts
