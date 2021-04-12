# coding: utf8


class Parse_clinical:
    """Creates the clinical files for the BIDS directory"""

    def __init__(self, path_clinical):
        import os.path as path

        import pandas as pd

        self.df_dict_mod = pd.read_csv(
            path.join(path_clinical, "clinical_info.tsv"), sep="\t"
        )
        self.df_clinical = pd.read_excel(
            path.join(path_clinical, "NIFD_Clinical_Data_2017_final_updated.xlsx")
        )
        self.df_ida = pd.read_csv(path.join(path_clinical, "ida.tsv"), sep="\t")

        self.merge_clinical = self.merge_clinical()
        self.df_clinical = self.merge_clinical

    def make_sessions_ida(self, pat_name):
        """Preprocesses ida for the left join operated in merge_clinical_scans.

        Args:
            pat_name: subject_ID
        """

        def write_ses(num):
            num = str(num)
            num = "0" + num if len(num) == 1 else num
            return "ses-M" + num

        bloc = self.df_ida[self.df_ida["Subject ID"] == pat_name]
        bloc = bloc[["Visit", "Study Date", "Age", "Weight", "Research Group"]]

        bloc["Visit"] = bloc["Visit"].apply(lambda x: int(x.split(" ")[1]))
        bloc = bloc.sort_values(["Visit"], ascending=[True])
        bloc["Visit"] = bloc["Visit"].apply(lambda x: write_ses(x))
        bloc = bloc.groupby("Visit").first().reset_index()
        bloc.columns = [
            "session_id",
            "examination_date",
            "age",
            "weight",
            "research_group",
        ]

        return bloc

    def merge_clinical_scans(self, pat_name):
        """Operate a left join between the ida and clinical table.

        For a given patient, we have the usual ida file, extended with the information from the clinical table.

        Args:
            pat_name: Name of the subject

        Returns:
            dfMerge: pandas dataframe corresponding to the left join
        """
        import warnings

        import pandas as pd

        bloc = self.df_clinical[self.df_clinical["LONI_ID"] == pat_name]

        def parse_date(x):
            x = str(x).split(" ")[0]
            x = x.split("-")
            sol = x[1] + "/" + x[2] + "/" + x[0]
            return sol

        warnings.simplefilter("ignore")
        bloc["CLINICAL_LINKDATE"] = bloc["CLINICAL_LINKDATE"].apply(
            lambda x: parse_date(x)
        )

        bloc_ida = self.make_sessions_ida(pat_name)
        bloc_ida["examination_date"] = bloc_ida["examination_date"].apply(
            lambda x: "0" + str(x).split(" ")[0]
            if str(x).split(" ")[0][0] != "0"
            else str(x).split(" ")[0]
        )

        dfMerge = pd.merge(
            bloc_ida, bloc, left_on="examination_date", right_on="CLINICAL_LINKDATE"
        )
        dfMerge = dfMerge.drop(columns=["CLINICAL_LINKDATE", "LONI_ID"])

        return dfMerge

    def merge_clinical(self):
        """Operate a left join between the clinical and ida table.

        The table that we obtain is the usual clinical table with data,
        when examination_date and CLINICAL_LINKDATE are the same, from the ida table.

        Returns:
            dfSol: pandas dataframe corresponding to the left join
        """
        # Could be optimized
        dfSol = self.df_clinical

        dfSol.insert(4, "Weight", "")
        dfSol.insert(4, "Research Group", "")
        dfSol.insert(4, "Age", "")

        def parse_date(x):
            x = str(x).split(" ")[0]
            x = x.split("-")
            sol = x[1] + "/" + x[2] + "/" + x[0]
            return sol

        curr_sub = None
        for index, row in dfSol.iterrows():
            if row["LONI_ID"] != curr_sub:
                curr_sub = row["LONI_ID"]
                dfMerge = self.merge_clinical_scans(curr_sub)
            info = dfMerge[
                dfMerge["examination_date"] == parse_date(row["CLINICAL_LINKDATE"])
            ]
            if not info.empty:
                dfSol.loc[index, "Weight"] = float(info["weight"])
                dfSol.loc[index, "Research Group"] = str(info.iloc[0]["research_group"])
                dfSol.loc[index, "Age"] = int(info["age"])

        return dfSol

    def get_clinical_ida(self):
        """Left join clinical data and ida, new version."""
        import math

        import pandas as pd

        df_clinical = self.df_clinical.copy()
        df_ida = self.df_ida.copy()
        df_ida["Session_number"] = df_ida["Visit"].apply(
            lambda x: int(x.split(" ")[-1])
        )

        df_ida = df_ida.groupby(["Subject ID", "Study Date"]).first().reset_index()
        df_ida = df_ida.sort_values(["Subject ID", "Session_number"]).reset_index()

        def parse_date(x):
            x = x.split("/")
            for i in range(len(x)):
                if len(x[i]) == 1:
                    x[i] = "0" + x[i]
            sol = x[1] + "/" + x[0] + "/" + x[2]
            return sol

        def parse_date2(x):
            x = str(x).split(" ")[0]
            x = x.split("-")
            sol = x[2] + "/" + x[1] + "/" + x[0]
            return sol

        df_ida["Study Date"] = df_ida["Study Date"].apply(lambda x: parse_date(x))
        df_clinical["CLINICAL_LINKDATE"] = df_clinical["CLINICAL_LINKDATE"].apply(
            lambda x: parse_date2(x)
        )

        dfSol = pd.merge(
            df_clinical,
            df_ida,
            how="left",
            left_on=["LONI_ID", "CLINICAL_LINKDATE"],
            right_on=["Subject ID", "Study Date"],
        )
        dfSol.insert(0, "session_id", "")

        dfSol["session_id"] = dfSol["Session_number"].apply(
            lambda x: ""
            if math.isnan(x)
            else ("ses-M" + str(int(x)) if x > 10 else "ses-M0" + str(int(x)))
        )

        return dfSol

    def make_sessions_type(self, pat_name, keep_all=False):
        """Updated version of make_sessions.

        Args:
            pat_name: subject_ID of a patient
            keep_all: if True, include clinical data not linked to a MRI, else include only clinical data linked to a MRI

        Returns:
            bloc: pandas dataframe corresponding to the "sessions.tsv" file
        """
        name_clinical, name_BIDS = self.get_names("sessions")
        name_clinical.insert(0, "session_id")
        name_BIDS.insert(0, "session_id")
        name_clinical.extend(["Age_x", "Research Group_x", "Weight_x"])
        name_BIDS.extend(["age", "research_group", "weight"])

        df_clinical_ida = self.get_clinical_ida()
        if not keep_all:
            df_clinical_ida = df_clinical_ida[df_clinical_ida["session_id"] != ""]

        bloc = df_clinical_ida[df_clinical_ida["LONI_ID"] == pat_name]
        bloc = bloc[name_clinical]
        bloc.columns = name_BIDS

        return bloc

    def get_names(self, file="sessions"):
        """
        Takes the clinical and BIDS_CLINICA fields that will be included in the clinical file "file"
        Returns them in list form

        Args:
            file: Corresponds to the clinical file that is about to be created, file = "sessions" or "participants"

        Returns:
            name_clinical, name_BIDS: 2 lists of fields
        """
        bloc = self.df_dict_mod[self.df_dict_mod["FILE"] == file]
        name_clinical = bloc["COLUMN_NAME"].tolist()
        name_BIDS = bloc["BIDS_CLINICA"].tolist()

        return name_clinical, name_BIDS

    def make_sessions(self, pat_name):
        """Create the sessions file for a given patient.

        Args:
            pat_name: subject_ID of a patient

        Returns:
            bloc: pandas dataframe corresponding to the "sessions.tsv" file
        """
        name_clinical, name_BIDS = self.get_names("sessions")

        bloc = self.df_clinical[self.df_clinical["LONI_ID"] == pat_name]
        bloc = bloc[name_clinical]
        bloc.columns = name_BIDS

        return bloc

    def make_participants(self, pat_list=None):
        """Create the participants file for all patients.

        Args:
            pat_list: subject_ID of all patients found in the converted BIDS directory

        Returns:
            bloc: pandas dataframe corresponding to the "participants.tsv" file
        """
        name_clinical, name_BIDS = self.get_names("participants")
        bloc = self.df_clinical.groupby("LONI_ID").first().reset_index()
        bloc = bloc[name_clinical]
        bloc["LONI_ID"] = bloc["LONI_ID"].apply(
            lambda x: "sub-NIFD" + x.replace("_", "")
        )
        bloc.columns = name_BIDS

        if pat_list is not None:
            bloc = bloc[bloc["participant_id"].isin(pat_list)]

        # Enforce the clinica_BIDS convention
        if "sex" in list(bloc):
            bloc["sex"] = bloc["sex"].apply(lambda x: "M" if x == 1 else "F")

        return bloc

    def make_scans(self, path_scans):
        """Create the scans file for a patient's session.

        Args:
            path_scans: path to a session for a patient

        Returns:
            bloc: pandas dataframe corresponding to the "participants.tsv" file
        """
        import os

        s = "filename	scan_id	mri_field\n"
        subs = [f.path.split("/")[-1] for f in os.scandir(path_scans) if f.is_dir()]
        for sub in subs:
            name = os.listdir(os.path.join(path_scans, sub))
            name = [
                i
                for i in name
                if i != ".DS_Store" and (i.endswith(".nii.gz") or i.endswith(".nii"))
            ]

            for n in name:
                s += sub + "/" + n + "\n"
        return s

    def write(self, df, path, name):
        """Save a pandas dataframe.

        Args:
            df: a pandas dataframe
            path: Path where the dataframe is to be saved
            name: name of the output file (Warning: do not include the extension, '.tsv' is added in the function)
        """
        import os

        name = os.path.join(path, name) + ".tsv"
        df = df.replace("", "n/a")
        df = df.fillna("n/a")
        df.to_csv(sep="\t", path_or_buf=name, index=False)

    def make_all(self, pathBIDS):
        """Make the participants.tsv and all sessions.tsv files for all subjects available in the BIDS directory.

        Args:
            pathBIDS: path to the BIDS directory
        """
        import os

        pat_list = os.listdir(pathBIDS)
        pat_list = [elt for elt in pat_list if elt.startswith("sub")]

        assert pat_list != [], "BIDS directory is empty"

        self.write(self.make_participants(pat_list), pathBIDS, "participants")

        for pat in pat_list:
            path_sessions = os.path.join(pathBIDS, pat)
            pat2 = pat[8] + "_S_" + pat[10:14]

            self.write(self.make_sessions_type(pat2), path_sessions, pat + "_sessions")

    def make_all_scans(self, to_convert):
        """Make the scans.tsv files for all subjects available in the BIDS directory.

        Args:
            to_convert: List of tuples of paths (path_in, path_out), computed for the initial image conversion
        """
        import os

        root = "/" + os.path.join(*to_convert[0][1].split("/")[:-4])

        def make_dic_tuples(to_convert, root):
            sol = {}
            dic_pat_sess = {}
            for tuple in to_convert:
                s_path1 = tuple[1].split("/")

                key = os.path.join(root, s_path1[-4], s_path1[-3])

                if s_path1[-4] not in dic_pat_sess:
                    dic_pat_sess[s_path1[-4]] = [s_path1[-3]]
                elif s_path1[-3] not in dic_pat_sess[s_path1[-4]]:
                    dic_pat_sess[s_path1[-4]].append(s_path1[-3])

                if key not in sol:
                    sol[key] = [tuple]
                else:
                    sol[key].append(tuple)

            return sol, dic_pat_sess

        def make_template():
            import pandas as pd

            new_cols = []

            for cell in list(self.df_ida["Imaging Protocol"]):
                if isinstance(cell, type("a")):
                    col_to_add = [i.split("=")[0] for i in cell.split(";")]
                    for i in col_to_add:
                        if i not in new_cols:
                            new_cols.append(i)

            new_cols.extend(["Modality", "Description", "Type", "Image ID"])
            new_cols.insert(0, "filename")
            sol = pd.DataFrame(columns=new_cols)
            return sol

        def extend_line(df_line_ida, template):
            col_values = {}
            s = list(df_line_ida["Imaging Protocol"])[0]
            for coup in s.split(";"):
                col_values[coup.split("=")[0]] = coup.split("=")[1]
            sol = df_line_ida
            for col_name in col_values:
                sol.insert(0, col_name, col_values[col_name])
            for name in list(template):
                if name not in list(sol):
                    sol.insert(0, name, "")
            return sol[list(template)]

        dic_tuples, dic_pat_sess = make_dic_tuples(to_convert, root)
        template = make_template()

        for sub_id in dic_pat_sess:
            for ses_num in dic_pat_sess[sub_id]:
                df_ses = template.copy()
                for tuple in dic_tuples[os.path.join(root, sub_id, ses_num)]:
                    s_path0 = tuple[0].split("/")
                    s_path1 = tuple[1].split("/")
                    filename = os.path.join(s_path1[-2], s_path1[-1]) + ".nii.gz"
                    df_line_ida = self.df_ida[
                        (
                            self.df_ida["Subject ID"]
                            == s_path1[-1][8] + "_S_" + s_path1[-1][10:14]
                        )
                        & (
                            self.df_ida["Visit"]
                            == "Month " + str(int(ses_num.split("M")[-1]))
                        )
                        & (self.df_ida["Description"] == s_path0[-3])
                    ]

                    if df_line_ida.empty:
                        df_line_ida = self.df_ida[
                            (
                                self.df_ida["Subject ID"]
                                == s_path1[-1][8] + "_S_" + s_path1[-1][10:14]
                            )
                            & (
                                self.df_ida["Visit"]
                                == "Month " + str(int(ses_num.split("M")[-1]))
                            )
                            & (
                                self.df_ida["Description"]
                                == s_path0[-3].replace("_", " ")
                            )
                        ]

                    # TR_BRAIN_3D_PIB_IR_CTAC -> TR:BRAIN 3D:PIB:IR CTAC
                    if df_line_ida.empty:
                        des = s_path0[-3].split("_")
                        if len(des) == 6:
                            des = (
                                f"{des[0]}:{des[1]} {des[2]}:{des[3]}:{des[4]} {des[5]}"
                            )
                        else:
                            des = f"{des[0]}:{des[1]} {des[2]}:{des[3]}:{des[4]} {des[5]} {des[6]}"
                        df_line_ida = self.df_ida[
                            (
                                self.df_ida["Subject ID"]
                                == s_path1[-1][8] + "_S_" + s_path1[-1][10:14]
                            )
                            & (
                                self.df_ida["Visit"]
                                == "Month " + str(int(ses_num.split("M")[-1]))
                            )
                            & (self.df_ida["Description"] == des)
                        ]

                    df_line_ida.insert(0, "filename", filename)
                    df_line_ida = extend_line(df_line_ida, template)

                    df_ses = df_ses.append(df_line_ida, ignore_index=True)

                columns = list(df_ses)
                columns[columns.index("Field Strength")] = "mri_field"
                df_ses.columns = columns

                self.write(
                    df_ses,
                    os.path.join(root, sub_id, ses_num),
                    sub_id + "_" + ses_num + "_scans",
                )
