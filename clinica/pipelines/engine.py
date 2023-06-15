"""This module contains the Pipeline abstract class needed for Clinica.

Subclasses are located in clinica/pipeline/<pipeline_name>/<pipeline_name>_pipeline.py
"""

import abc
from typing import Tuple

import click
from nipype.pipeline.engine import Workflow


def clinica_pipeline(func):
    """Turns a CLI implementation into a Clinica Pipeline.

    More precisely, this decorator takes care of:
        - registering the implemented command as a subcommand
          of `clinica run`.
    """
    from clinica.pipelines.cli import cli as run_cli

    run_cli.add_command(func)

    return func


def postset(attribute, value):
    """Sets the attribute of an object after the execution.

    Args:
        attribute: An object's attribute to be set.
        value: A desired value for the object's attribute.

    Returns:
        A decorator executed after the decorated function is.
    """

    def postset_decorator(func):
        def func_wrapper(self, *args, **kwargs):
            res = func(self, *args, **kwargs)
            setattr(self, attribute, value)
            return res

        return func_wrapper

    return postset_decorator


def detect_cross_sectional_and_longitudinal_subjects(
    all_subs: list, bids_dir: str
) -> Tuple[list, list]:
    """This function detects whether a subject has a longitudinal or cross sectional structure.
    If it has a mixed structure, it will be considered cross sectionnal.

    Parameters
    ----------
    all_subs: list
        List containing all the names of the subjects.

    bids_dir: str
        Path to the bids directory.

    Returns
    -------
    cross_subj: list
        List of all the subject detected with a cross_sectional organisation.

    long_subj: list
        List of all the subject detected with a longitudinal organisation.

    Examples
    --------
    BIDS
    ├── sub-01
    │   └── anat
    └── sub-02
        └── ses-M000
    >>> detect_cross_sectional_and_longitudinal_subjects(["sub-01", "sub-02"], /Users/name.surname/BIDS)
    (["sub-01"], ["sub-02])
    """
    from os import listdir
    from os.path import isdir, join

    cross_subj = []
    long_subj = []
    for sub in all_subs:
        folder_list = [
            f for f in listdir(join(bids_dir, sub)) if isdir(join(bids_dir, sub, f))
        ]

        if not all([fold.startswith("ses-") for fold in folder_list]):
            cross_subj.append(sub)
        else:
            long_subj.append(sub)
    return cross_subj, long_subj


class Pipeline(Workflow):
    """Clinica Pipeline class.

    This class overwrites the `Workflow` to integrate and encourage the
    use of BIDS and CAPS data structures as inputs and outputs of the pipelines
    developed for the Clinica software.

    The global architecture of a Clinica pipelines is as follow:
        [ Data Input Stream ]
                |
            [ Input ]
                |
            [[ Core ]] <- Could be one or more `npe.Node`s
                |
            [ Output ]
                |
        [ Data Output Stream ]

    Attributes:
        is_built (bool): Informs if the pipelines has been built or not.
        parameters (dict): Parameters of the pipelines.
        info (dict): Information presented in the associated `info.json` file.
        input_node (:obj:`npe.Node`): Identity interface connecting inputs.
        output_node (:obj:`npe.Node`): Identity interface connecting outputs.
        bids_directory (str): Directory used to read the data from, in BIDS.
        caps_directory (str): Directory used to read/write the data from/to,
            in CAPS.
        subjects (list): List of subjects defined in the `subjects.tsv` file.
            # TODO(@jguillon): Check the subjects-sessions file name.
        sessions (list): List of sessions defined in the `subjects.tsv` file.
        tsv_file (str): Path to the subjects-sessions `.tsv` file.
        info_file (str): Path to the associated `info.json` file.
    """

    __metaclass__ = abc.ABCMeta

    def __init__(
        self,
        bids_directory=None,
        caps_directory=None,
        tsv_file=None,
        overwrite_caps=False,
        base_dir=None,
        parameters={},
        name=None,
    ):
        """Init a Pipeline object.

        Args:
            bids_directory (str, optional): Path to a BIDS directory. Defaults to None.
            caps_directory (str, optional): Path to a CAPS directory. Defaults to None.
            tsv_file (str, optional): Path to a subjects-sessions `.tsv` file. Defaults to None.
            overwrite_caps (bool, optional): Overwrite or not output directory.. Defaults to False.
            base_dir (str, optional): Working directory (attribute of Nipype::Workflow class). Defaults to None.
            parameters (dict, optional): Pipeline parameters. Defaults to {}.
            name (str, optional): Pipeline name. Defaults to None.

        Raises:
            RuntimeError: [description]
        """
        import inspect
        import os
        from tempfile import mkdtemp

        from clinica.utils.inputs import check_bids_folder, check_caps_folder
        from clinica.utils.participant import get_subject_session_list

        self._is_built = False
        self._overwrite_caps = overwrite_caps
        self._bids_directory = bids_directory
        self._caps_directory = caps_directory
        self._verbosity = "debug"
        self._tsv_file = tsv_file
        self._info_file = os.path.join(
            os.path.dirname(os.path.abspath(inspect.getfile(self.__class__))),
            "info.json",
        )
        self._info = {}

        if base_dir:
            self.base_dir = base_dir
            self._base_dir_was_specified = True
        else:
            self.base_dir = mkdtemp()
            self._base_dir_was_specified = False

        self._name = name or self.__class__.__name__
        self._parameters = parameters

        if not self._bids_directory:
            if not self._caps_directory:
                raise RuntimeError(
                    f"The {self._name} pipeline does not contain "
                    "BIDS nor CAPS directory at the initialization."
                )

            check_caps_folder(self._caps_directory)
            input_dir = self._caps_directory
            is_bids_dir = False
        else:
            check_bids_folder(self._bids_directory)
            input_dir = self._bids_directory
            is_bids_dir = True
        self._sessions, self._subjects = get_subject_session_list(
            input_dir, tsv_file, is_bids_dir, False, base_dir
        )

        self.init_nodes()

    @staticmethod
    def get_processed_images(caps_directory, subjects, sessions):
        """Extract processed image IDs in `caps_directory` based on `subjects`_`sessions`.

        Todo:
            [ ] Implement this static method in all pipelines
            [ ] Make it abstract to force overload in future pipelines
        """
        from clinica.utils.exceptions import ClinicaException
        from clinica.utils.stream import cprint

        cprint(msg="Pipeline finished with errors.", lvl="error")
        cprint(msg="CAPS outputs were not found for some image(s):", lvl="error")
        raise ClinicaException(
            "Implementation on which image(s) failed will appear soon."
        )

    def init_nodes(self):
        """Init the basic workflow and I/O nodes necessary before build."""
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        if self.get_input_fields():
            self._input_node = npe.Node(
                name="Input",
                interface=nutil.IdentityInterface(
                    fields=self.get_input_fields(), mandatory_inputs=False
                ),
            )
        else:
            self._input_node = None

        if self.get_output_fields():
            self._output_node = npe.Node(
                name="Output",
                interface=nutil.IdentityInterface(
                    fields=self.get_output_fields(), mandatory_inputs=False
                ),
            )
        else:
            self._output_node = None

        Workflow.__init__(self, self._name, self.base_dir)
        if self.input_node:
            self.add_nodes([self.input_node])
        if self.output_node:
            self.add_nodes([self.output_node])

    def has_input_connections(self):
        """Checks if the Pipeline's input node has been connected.

        Returns:
            True if the input node is connected, False otherwise.
        """
        if self.input_node:
            return self._graph.in_degree(self.input_node) > 0
        else:
            return False

    def has_output_connections(self):
        """Checks if the Pipeline's output node has been connected.

        Returns:
            True if the output node is connected, False otherwise.
        """
        if self.output_node:
            return self._graph.out_degree(self.output_node) > 0
        else:
            return False

    @postset("is_built", True)
    def build(self):
        """Builds the core, input and output nodes of the Pipeline.

        This method first checks it has already been run. It then checks
        the pipelines dependencies and, in this order, builds the core nodes,
        the input node and, finally, the output node of the Pipeline.

        Since this method returns the concerned object, it can be chained to
        any other method of the Pipeline class.

        Returns:
            self: A Pipeline object.
        """
        if not self.is_built:
            self.check_dependencies()
            self.check_pipeline_parameters()
            if not self.has_input_connections():
                self.build_input_node()
            self.build_core_nodes()
            if not self.has_output_connections():
                self.build_output_node()
        return self

    def run(self, plugin=None, plugin_args=None, update_hash=False, bypass_check=False):
        """Executes the Pipeline.

        It overwrites the default Workflow method to check if the
        Pipeline is built before running it. If not, it builds it and then
        run it.
        It also checks whether there is enough space left on the disks, and if
        the number of threads to run in parallel is consistent with what is
        possible on the CPU.

        Args:
            Similar to those of Workflow.run.

        Returns:
            An execution graph (see Workflow.run).
        """
        import shutil

        from networkx import Graph, NetworkXError

        from clinica.utils.stream import cprint
        from clinica.utils.ux import print_failed_images

        if not self.is_built:
            self.build()
        self.check_not_cross_sectional()

        if not bypass_check:
            self.check_size()
            plugin_args = self.update_parallelize_info(plugin_args)
            plugin = "MultiProc"
        exec_graph = []
        try:
            exec_graph = Workflow.run(self, plugin, plugin_args, update_hash)
            if not self.base_dir_was_specified:
                shutil.rmtree(self.base_dir)

        except RuntimeError as e:
            # Check that it is a Nipype error
            if "Workflow did not execute cleanly. Check log for details" in str(e):
                input_ids = [
                    p_id + "_" + s_id
                    for p_id, s_id in zip(self.subjects, self.sessions)
                ]
                output_ids = self.get_processed_images(
                    caps_directory=self.caps_directory,
                    subjects=self.subjects,
                    sessions=self.sessions,
                )
                missing_ids = list(set(input_ids) - set(output_ids))
                print_failed_images(self.name, missing_ids)
            else:
                raise e
        except NetworkXError:
            cprint(
                msg=(
                    "Either all the images were already run by the pipeline "
                    "or no image was found to run the pipeline."
                ),
                lvl="warning",
            )
            exec_graph = Graph()
        return exec_graph

    def load_info(self):
        """Loads the associated info.json file.

        Todo:
            - [ ] Raise an appropriate exception when the info file can't open

        Raises:
            None. # TODO(@jguillon)

        Returns:
            self: A Pipeline object.
        """
        import json

        with open(self.info_file) as info_file:
            self.info = json.load(info_file)
        return self

    def check_dependencies(self):
        """Checks if listed dependencies are present.

        Loads the pipelines related `info.json` file and check each one of the
        dependencies listed in the JSON "dependencies" field. Its raises
        exception if a program in the list does not exist or if environment
        variables are not properly defined.

        Todo:
            - [ ] MATLAB toolbox dependency checking
            - [x] check MATLAB
            - [ ] Clinica pipelines dependency checks
            - [ ] Check dependencies version

        Raises:
            Exception: Raises an exception when bad dependency types given in
            the `info.json` file are detected.

        Returns:
            self: A Pipeline object.
        """
        import clinica.utils.check_dependency as chk

        # Checking functions preparation
        check_software = {
            # 'matlab': chk.check_matlab,
            "ants": chk.check_ants,
            "convert3d": chk.check_convert3d,
            "spm": chk.check_spm,
            "freesurfer": chk.check_freesurfer,
            "fsl": chk.check_fsl,
            "mrtrix": chk.check_mrtrix,
            "matlab": chk.check_matlab,
            "petpvc": chk.check_petpvc,
        }
        check_binary = chk.is_binary_present
        # check_toolbox = chk.is_toolbox_present
        # check_pipeline = chk.is_pipeline_present

        # Load the info.json file
        if not self.info:
            self.load_info()

        # Dependencies checking
        for d in self.info["dependencies"]:
            if d["type"] == "software":
                check_software[d["name"]]()
            elif d["type"] == "binary":
                check_binary(d["name"])
            elif d["type"] == "toolbox":
                pass
            elif d["type"] == "pipeline":
                pass
            else:
                raise Exception(
                    f"Pipeline.check_dependencies() Unknown dependency type: '{d['type']}'."
                )

        self.check_custom_dependencies()

        return self

    def check_size(self):
        """Check if the pipeline has enough space on the disk for both working directory and CAPS."""
        import select
        import sys
        from os import statvfs
        from os.path import dirname

        from clinica.utils.stream import cprint

        SYMBOLS = {
            "customary": ("B", "K", "M", "G", "T", "P", "E", "Z", "Y"),
            "customary_ext": (
                "byte",
                "kilo",
                "mega",
                "giga",
                "tera",
                "peta",
                "exa",
                "zetta",
                "iotta",
            ),
            "iec": ("Bi", "Ki", "Mi", "Gi", "Ti", "Pi", "Ei", "Zi", "Yi"),
            "iec_ext": (
                "byte",
                "kibi",
                "mebi",
                "gibi",
                "tebi",
                "pebi",
                "exbi",
                "zebi",
                "yobi",
            ),
        }

        def bytes2human(n, format="%(value).1f %(symbol)s", symbols="customary"):
            """
            Convert n bytes into a human readable string based on format.
            Symbols can be either "customary", "customary_ext", "iec" or "iec_ext",
            see: http://goo.gl/kTQMs

            License:
            Bytes-to-human / human-to-bytes converter.
            Based on: http://goo.gl/kTQMs
            Working with Python 2.x and 3.x.
            Author: Giampaolo Rodola' <g.rodola [AT] gmail [DOT] com>
            License: MIT
            """
            n = int(n)
            if n < 0:
                raise ValueError("n < 0")
            symbols = SYMBOLS[symbols]
            prefix = {}
            for i, s in enumerate(symbols[1:]):
                prefix[s] = 1 << (i + 1) * 10
            for symbol in reversed(symbols[1:]):
                if n >= prefix[symbol]:
                    value = float(n) / prefix[symbol]
                    return format % locals()
            return format % dict(symbol=symbols[0], value=n)

        def human2bytes(s):
            """
            Attempts to guess the string format based on default symbols
            set and return the corresponding bytes as an integer.
            When unable to recognize the format ValueError is raised.

            License:
            Bytes-to-human / human-to-bytes converter.
            Based on: http://goo.gl/kTQMs
            Working with Python 2.x and 3.x.
            Author: Giampaolo Rodola' <g.rodola [AT] gmail [DOT] com>
            License: MIT
            """
            init = s
            num = ""
            while s and s[0:1].isdigit() or s[0:1] == ".":
                num += s[0]
                s = s[1:]
            num = float(num)
            letter = s.strip()
            for name, sset in SYMBOLS.items():
                if letter in sset:
                    break
            else:
                if letter == "k":
                    # treat 'k' as an alias for 'K' as per: http://goo.gl/kTQMs
                    sset = SYMBOLS["customary"]
                    letter = letter.upper()
                else:
                    raise ValueError("can't interpret %r" % init)
            prefix = {sset[0]: 1}
            for i, s in enumerate(sset[1:]):
                prefix[s] = 1 << (i + 1) * 10
            return int(num * prefix[letter])

        timeout = 15

        # Get the number of sessions
        n_sessions = len(self.subjects)
        try:
            caps_stat = statvfs(self.caps_directory)
        except FileNotFoundError:
            # CAPS folder may not exist yet
            caps_stat = statvfs(dirname(self.caps_directory))
        wd_stat = statvfs(dirname(self.base_dir))

        # Estimate space left on partition/disk/whatever caps and wd is located
        free_space_caps = caps_stat.f_bavail * caps_stat.f_frsize
        free_space_wd = wd_stat.f_bavail * wd_stat.f_frsize

        try:
            space_needed_caps_1_session = self.info["space_caps"]
            space_needed_wd_1_session = self.info["space_wd"]
            space_needed_caps = n_sessions * human2bytes(space_needed_caps_1_session)
            space_needed_wd = n_sessions * human2bytes(space_needed_wd_1_session)
            error = ""
            if free_space_caps == free_space_wd:
                if space_needed_caps + space_needed_wd > free_space_wd:
                    # We assume this is the same disk
                    error = (
                        f"{error}Space needed for CAPS and working directory ("
                        f"{bytes2human(space_needed_caps + space_needed_wd)}"
                        f") is greater than what is left on your hard drive ("
                        f"{bytes2human(free_space_wd)})"
                    )
            else:
                if space_needed_caps > free_space_caps:
                    error = (
                        f"{error} Space needed for CAPS ("
                        f"{bytes2human(space_needed_caps)}"
                        f") is greater than what is left on your hard drive ("
                        f"{bytes2human(free_space_caps)})\n"
                    )
                if space_needed_wd > free_space_wd:
                    error = (
                        f"{error}Space needed for working_directory ("
                        f"{bytes2human(space_needed_wd)}"
                        f") is greater than what is left on your hard drive ("
                        f"{bytes2human(free_space_wd)})\n"
                    )
            if error:
                cprint(msg=f"[SpaceError] {error}", lvl="error")
                if not click.confirm(
                    "Space issues detected. Do you still want to run the pipeline?"
                ):
                    click.echo("Clinica will now exit...")
                    sys.exit()

        except KeyError:
            cprint("No info on how much size the pipeline takes. Running anyway...")

    def update_parallelize_info(self, plugin_args):
        """Performs some checks of the number of threads given in parameters,
        given the number of CPUs of the machine in which clinica is running.
        We force the use of plugin MultiProc

        Author: Arnaud Marcoux"""
        import select
        import sys
        from multiprocessing import cpu_count

        from clinica.utils.stream import cprint

        # count number of CPUs
        n_cpu = cpu_count()
        # timeout value: max time allowed to decide how many thread
        # to run in parallel (sec)
        timeout = 15

        # Use this var to know in the end if we need to ask the user
        # an other number
        ask_user = False

        try:
            # if no --n_procs arg is used, plugin_arg is None
            # so we need a try / except block
            n_thread_cmdline = plugin_args["n_procs"]
            if n_thread_cmdline > n_cpu:
                cprint(
                    msg=(
                        f"You are trying to run clinica with a number of threads ({n_thread_cmdline}) superior to your "
                        f"number of CPUs ({n_cpu})."
                    ),
                    lvl="warning",
                )
                ask_user = True
        except TypeError:
            cprint(
                msg=f"You did not specify the number of threads to run in parallel (--n_procs argument).",
                lvl="warning",
            )
            cprint(
                msg=(
                    f"Computation time can be shorten as you have {n_cpu} CPUs on this computer. "
                    f"We recommend using {n_cpu-1} threads."
                ),
                lvl="warning",
            )
            ask_user = True

        if ask_user:
            n_procs = click.prompt(
                text="How many threads do you want to use?",
                default=max(1, n_cpu - 1),
                show_default=True,
            )
            click.echo(
                f"Number of threads set to {n_procs}. You may set the --n_procs argument "
                "to disable this message for future calls."
            )
            if plugin_args:
                plugin_args["n_procs"] = n_procs
            else:
                plugin_args = {"n_procs": n_procs}

        return plugin_args

    def check_not_cross_sectional(self):
        """
        This function checks if the dataset is longitudinal. If it is cross
        sectional, clinica proposes to convert it in a clinica compliant form.

        author: Arnaud Marcoux
        """
        import sys
        from os import listdir
        from os.path import abspath, basename, dirname, isdir, join

        from clinica.utils.exceptions import ClinicaInconsistentDatasetError
        from clinica.utils.stream import cprint

        def _check_cross_subj(cross_subj: list) -> None:
            if len(cross_subj) > 0:
                raise ClinicaInconsistentDatasetError(cross_subj)

        def convert_cross_sectional(bids_in, bids_out, cross_subjects, long_subjects):
            """
            This function converts a cross-sectional-bids dataset into a
            longitudinal clinica-compliant dataset

            Args:
                bids_in: cross sectional bids dataset you want to convert
                bids_out: converted longitudinal bids dataset
                cross_subjects: list of subjects in cross sectional form
                (they need some adjustment)
                long_subjects: list of subjects in longitudinal form (they
                just need to be copied)

            Returns:
                nothing
            """
            from os import mkdir
            from os.path import exists, isfile
            from shutil import copy2, copytree

            def add_ses(f):
                """
                Use this function to transform a cross sectional filename into
                a longitudinal one.
                Examples:
                sub-ADNI001_scans.tsv -> sub-ADNI001_ses-M000_scans.tsv
                sub-023a_ses-M012_T1w.nii.gz -> sub-023a_ses-M012_T1w.nii.gz (no
                    modification done if filename already has a session)

                Args:
                    f: filename

                Returns:
                    filename with '_ses-M000_ added just after participant_id
                """
                import re

                # If filename contains ses-..., returns the original filename
                # Regex explication:
                # ^ start of string
                # ([a-zA-Z0-9]*) matches any number of characters from a to z,
                #       A to Z, 0 to 9, and store it in group(1)
                # (?!ses-[a-zA-Z0-9]) do not match if there is already a 'ses-'
                # (.*) catches the rest of the string
                m = re.search(r"(^sub-[a-zA-Z0-9]*)_(?!ses-[a-zA-Z0-9])(.*)", f)
                try:
                    return m.group(1) + "_ses-M000_" + m.group(2)
                except AttributeError:
                    # If something goes wrong, we return the original filename
                    return f

            def copy2_add_ses(src, dst):
                """
                copy2_add_ses calls copy2 function from shutil, but modifies
                the filename of the copied files if they match the regex
                template described in add_ses() function

                Args:
                    src: path to the file that needs to be copied
                    dst: original destination for the copied file

                Returns:
                    copy2 with modified filename
                """
                from os.path import basename, dirname, join
                from shutil import copy2

                dst_modified = join(dirname(dst), add_ses(basename(src)))
                return copy2(src, dst_modified)

            if not exists(bids_out):
                # Create the output folder if it does not exists yet
                mkdir(bids_out)

            # First part of the algorithm: deal with subjects that does not
            # have longitudinal (session) information
            for subj in cross_subjects:
                # Get list of of files/folders to copy. Remove hidden element
                # ( if they start with a dot '.')
                to_copy = [
                    f for f in listdir(join(bids_in, subj)) if not f.startswith(".")
                ]
                # Always check that folder are existing, otherwise an exception
                # is raised even though that does not prevent the function from
                # working
                if not exists(join(bids_out, subj)):
                    mkdir(join(bids_out, subj))
                if not exists(join(bids_out, subj, "ses-M000")):
                    mkdir(join(bids_out, subj, "ses-M000"))
                for el in to_copy:
                    path_el = join(bids_in, subj, el)
                    if not exists(join(bids_out, subj, "ses-M000", el)):
                        # If the element to copy is a folder...
                        if isdir(path_el):
                            # Copytree is used with a 'custom' copy function,
                            # to give a correct filemename to the files inside
                            # the copied folders
                            copytree(
                                path_el,
                                join(bids_out, subj, "ses-M000", basename(path_el)),
                                copy_function=copy2_add_ses,
                            )
                        # If the element to copy is a file...
                        elif isfile(path_el):
                            # Modify the filename with add_ses function
                            new_filename_wo_ses = add_ses(el)
                            copy2(
                                path_el,
                                join(bids_out, subj, "ses-M000", new_filename_wo_ses),
                            )
            # Second part of the algorithm: deal with subjects that do not
            # have the problem. We only xopy the content of the folder, and no
            # filename needs to be changed
            for su in long_subjects:
                # Do not copy hidden files
                to_copy = [
                    f for f in listdir(join(bids_in, su)) if not f.startswith(".")
                ]
                if not exists(join(bids_out, su)):
                    mkdir(join(bids_out, su))
                for el in to_copy:
                    path_el = join(bids_in, su, el)
                    if not exists(join(bids_out, su, el)):
                        # 2 possible cases: element to pcopy is a folder or a
                        # file
                        if isdir(path_el):
                            copytree(path_el, join(bids_out, su, basename(path_el)))
                        elif isfile(path_el):
                            copy2(path_el, join(bids_out, su))

        if self.bids_directory is not None:
            bids_dir = abspath(self.bids_directory)
            # Extract all subjects in BIDS directory: element must be a folder and
            # a name starting with 'sub-'
            all_subs = [
                f
                for f in listdir(bids_dir)
                if isdir(join(bids_dir, f)) and f.startswith("sub-")
            ]

            cross_subj, long_subj = detect_cross_sectional_and_longitudinal_subjects(
                all_subs, bids_dir
            )
            try:
                _check_cross_subj(cross_subj)
            except ClinicaInconsistentDatasetError as e:
                cprint(e, lvl="warning")
                proposed_bids = join(
                    dirname(bids_dir), basename(bids_dir) + "_clinica_compliant"
                )

                if not click.confirm(
                    "Do you want to proceed with the conversion in another "
                    "folder? (Your original BIDS folder will not be modified "
                    f"and the folder {proposed_bids} will be created.)"
                ):
                    click.echo("Clinica will now exit...")
                    sys.exit()
                else:
                    cprint("Converting cross-sectional dataset into longitudinal...")
                    convert_cross_sectional(
                        bids_dir, proposed_bids, cross_subj, long_subj
                    )
                    cprint(
                        f"Conversion succeeded. Your clinica-compliant dataset is located here: {proposed_bids}"
                    )

    @property
    def base_dir_was_specified(self):
        return self._base_dir_was_specified

    @property
    def is_built(self):
        return self._is_built

    @is_built.setter
    def is_built(self, value):
        self._is_built = value

    @property
    def overwrite_caps(self):
        return self._overwrite_caps

    @property
    def parameters(self):
        return self._parameters

    @parameters.setter
    def parameters(self, value):
        self._parameters = value
        # Need to rebuild input, output and core nodes
        self.is_built = False
        self.init_nodes()

    @property
    def info(self):
        return self._info

    @info.setter
    def info(self, value):
        self._info = value

    @property
    def input_node(self):
        return self._input_node

    @property
    def output_node(self):
        return self._output_node

    @property
    def bids_directory(self):
        return self._bids_directory

    @property
    def caps_directory(self):
        return self._caps_directory

    @property
    def subjects(self):
        return self._subjects

    @subjects.setter
    def subjects(self, value):
        self._subjects = value
        self.is_built = False

    @property
    def sessions(self):
        return self._sessions

    @sessions.setter
    def sessions(self, value):
        self._sessions = value
        self.is_built = False

    @property
    def tsv_file(self):
        return self._tsv_file

    @property
    def info_file(self):
        return self._info_file

    @abc.abstractmethod
    def build_core_nodes(self):
        """Builds the Pipeline's core nodes.

        This method should use and connect to the `Pipeline.input_node` one
        or more core `Node`s. The outputs of the core processing should then be
        connected to the `Pipeline.output_node`.
        """

    @abc.abstractmethod
    def build_input_node(self):
        """Builds the Pipeline's input data stream node.

        Warnings:
            This method does not modify the `Pipeline.input_node` (see the
            notes about the global architecture in the class documentation).
        """

    @abc.abstractmethod
    def build_output_node(self):
        """Builds the Pipeline's output data stream node.

        Warnings:
            This method does not modify the `Pipeline.output_node` (see the
            notes about the global architecture in the class documentation).
        """

    @abc.abstractmethod
    def get_input_fields(self):
        """Lists the input fields of the Pipeline.

        Returns:
            A list of strings defining the fields of the `IdentityInterface`
            of the `Pipeline.input_node`.
        """

    @abc.abstractmethod
    def get_output_fields(self):
        """Lists the output fields of the Pipeline.

        Returns:
            A list of strings defining the fields of the `IdentityInterface`
            of the `Pipeline.output_node`.
        """

    @abc.abstractmethod
    def check_custom_dependencies(self):
        """Check dependencies provided by the developer."""

    @abc.abstractmethod
    def check_pipeline_parameters(self):
        """Check pipeline parameters."""


class PETPipeline(Pipeline):
    def check_pipeline_parameters(self):
        """Check pipeline parameters."""
        if "acq_label" not in self.parameters.keys():
            raise KeyError("Missing compulsory acq_label key in pipeline parameter.")
        self.parameters.setdefault("reconstruction_method", None)

    def _get_pet_scans_query(self) -> dict:
        """Return the query to retrieve PET scans."""
        from clinica.utils.input_files import bids_pet_nii
        from clinica.utils.pet import ReconstructionMethod, Tracer

        pet_tracer = None
        if self.parameters["acq_label"] is not None:
            pet_tracer = Tracer(self.parameters["acq_label"])

        reconstruction_method = None
        if self.parameters["reconstruction_method"] is not None:
            reconstruction_method = ReconstructionMethod(
                self.parameters["reconstruction_method"]
            )

        return bids_pet_nii(pet_tracer, reconstruction_method)
