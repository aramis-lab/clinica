"""This module contains the Pipeline abstract class needed for Clinica.

Subclasses are located in clinica/pipeline/<pipeline_name>/<pipeline_name>_pipeline.py
"""

import abc
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

import click
from nipype.interfaces.utility import IdentityInterface
from nipype.pipeline.engine import Node, Workflow

from clinica.utils.bids import Visit
from clinica.utils.check_dependency import SoftwareDependency, ThirdPartySoftware
from clinica.utils.group import GroupID, GroupLabel
from clinica.utils.stream import log_and_warn


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


def _detect_cross_sectional_and_longitudinal_subjects(
    subjects: Iterable[str], bids_dir: Path
) -> Tuple[List[str], List[str]]:
    """Detect whether a subject's folder has a longitudinal or cross-sectional structure.

    If it has a mixed structure, it will be considered cross-sectionnal.

    Parameters
    ----------
    subjects: Iterable of str
        List containing all the names of the subjects.

    bids_dir: Path
        The path to the bids directory.

    Returns
    -------
    cross_sectional_subjects: list of str
        List of all the subject detected with a cross_sectional organisation.

    longitudinal_subjects: list of str
        List of all the subject detected with a longitudinal organisation.

    Examples
    --------
    BIDS
    ├── sub-01
    │   └── anat
    └── sub-02
        └── ses-M000
    >>> _detect_cross_sectional_and_longitudinal_subjects(["sub-01", "sub-02"], Path("/Users/name.surname/BIDS"))
    (["sub-01"], ["sub-02])
    """
    cross_sectional_subjects = []
    longitudinal_subjects = []
    for subject in subjects:
        folders = [f.name for f in (bids_dir / subject).iterdir() if f.is_dir()]
        if not all([folder.startswith("ses-") for folder in folders]):
            cross_sectional_subjects.append(subject)
        else:
            longitudinal_subjects.append(subject)
    return cross_sectional_subjects, longitudinal_subjects


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


def _bytes2human(
    n: int, format: str = "%(value).1f %(symbol)s", symbols: str = "customary"
) -> str:
    """Convert n bytes into a human readable string based on format.

    Symbols can be either "customary", "customary_ext", "iec" or "iec_ext", see: http://goo.gl/kTQMs

    License:
    Bytes-to-human / human-to-bytes converter.
    Based on: http://goo.gl/kTQMs
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


def _human2bytes(s: str) -> int:
    """Attempts to guess the string format based on default symbols
    set and return the corresponding bytes as an integer.

    When unable to recognize the format ValueError is raised.

    License:
    Bytes-to-human / human-to-bytes converter.
    Based on: http://goo.gl/kTQMs
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


def _add_session_label(filename: str) -> str:
    """Transform a cross-sectional filename into a longitudinal one.

    Examples
    --------
    >>> _add_session_label("sub-ADNI001_scans.tsv")
    sub-ADNI001_ses-M000_scans.tsv

    No modification done if filename already has a session:

    >>> _add_session_label("sub-023a_ses-M012_T1w.nii.gz")
    sub-023a_ses-M012_T1w.nii.gz

    Parameters
    ----------
    filename : str

    Returns
    -------
    str :
        filename with '_ses-M000_ added just after participant_id
    """
    import re

    m = re.search(r"(^sub-[a-zA-Z0-9]*)_(?!ses-[a-zA-Z0-9])(.*)", filename)
    try:
        return m.group(1) + "_ses-M000_" + m.group(2)
    except AttributeError:
        return filename


def _copy_and_add_session_label(src: str, dst: str) -> str:
    """Copy src to dst, but modifies the filename of the copied files if they
    match the regex template described in the _add_session_label() function.

    Parameters
    ----------
    src : str
        The path to the file that needs to be copied.

    dst : str
        The original destination for the copied file.

    Returns
    -------
    str :
        Copy with modified filename.
    """
    from shutil import copy2

    filename = Path(src).name

    return copy2(src, str(Path(dst).parent / _add_session_label(filename)))


def _convert_cross_sectional(
    bids_in: Path,
    bids_out: Path,
    cross_subjects: Iterable[str],
    long_subjects: Iterable[str],
):
    """Convert a cross-sectional-bids dataset into a longitudinal clinica-compliant dataset.

    Parameters
    ----------
    bids_in : Path
        The path to the cross-sectional bids dataset you want to convert.

    bids_out : Path
        The path to the converted longitudinal bids dataset.

    cross_subjects : Iterable of str
        Subjects in cross-sectional form (they need some adjustment).

    long_subjects : Iterable of str
        Subjects in longitudinal form (they just need to be copied).
    """
    if not bids_out.exists():
        bids_out.mkdir()
    _copy_metadata_files(bids_in, bids_out, ("dataset_description.json", ".bidsignore"))
    _convert(bids_in, bids_out, cross_subjects, cross_sectional=True)
    _convert(bids_in, bids_out, long_subjects, cross_sectional=False)


def _copy_metadata_files(bids_in: Path, bids_out: Path, files_to_copy: Iterable[str]):
    from shutil import copy2

    for file in files_to_copy:
        if (bids_in / file).is_file():
            copy2(bids_in / file, bids_out / file)


def _convert(
    bids_in: Path, bids_out: Path, subjects: Iterable[str], cross_sectional: bool
):
    for subject in subjects:
        output_folder = (
            bids_out / subject / "ses-M000" if cross_sectional else bids_out / subject
        )
        if not output_folder.exists():
            output_folder.mkdir(parents=True)
        for file_to_copy in [
            f for f in (bids_in / subject).iterdir() if not f.name.startswith(".")
        ]:
            copy_func = _copy_cross_sectional if cross_sectional else _copy_longitudinal
            copy_func(bids_in / subject / file_to_copy, output_folder)


def _copy_cross_sectional(file_to_copy: Path, output_folder: Path):
    from shutil import copy2, copytree

    if not (output_folder / file_to_copy.name).exists():
        if file_to_copy.is_dir():
            copytree(
                file_to_copy,
                output_folder / file_to_copy.name,
                copy_function=_copy_and_add_session_label,
            )
        elif file_to_copy.is_file():
            new_filename_wo_ses = _add_session_label(file_to_copy.name)
            copy2(file_to_copy, output_folder / new_filename_wo_ses)


def _copy_longitudinal(file_to_copy: Path, output_folder: Path) -> None:
    from shutil import copy2, copytree

    if not (output_folder / file_to_copy.name).exists():
        if file_to_copy.is_dir():
            copytree(file_to_copy, output_folder / file_to_copy.name)
        elif file_to_copy.is_file():
            copy2(file_to_copy, output_folder)


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

    Attributes
    ----------
    is_built : bool
        Informs if the pipelines has been built or not.

    parameters : dict
        Parameters of the pipelines.

    info : dict
        Information presented in the associated `info.json` file.

    input_node : :obj:`npe.Node`
        Identity interface connecting inputs.

    output_node : :obj:`npe.Node`
        Identity interface connecting outputs.

    bids_directory : Path
        Directory used to read the data from, in BIDS.

    caps_directory : Path
        Directory used to read/write the data from/to in CAPS.

    subjects : list
        List of subjects defined in the `subjects.tsv` file.
        TODO(@jguillon): Check the subjects-sessions file name.

    sessions : list
        List of sessions defined in the `subjects.tsv` file.

    tsv_file : str
        Path to the subjects-sessions `.tsv` file.

    info_file : str
        Path to the associated `info.json` file.
    """

    __metaclass__ = abc.ABCMeta

    def __init__(
        self,
        bids_directory: Optional[str] = None,
        caps_directory: Optional[str] = None,
        tsv_file: Optional[str] = None,
        overwrite_caps: Optional[bool] = False,
        base_dir: Optional[str] = None,
        parameters: Optional[dict] = None,
        name: Optional[str] = None,
        ignore_dependencies: Optional[List[str]] = None,
        caps_name: Optional[str] = None,
    ):
        """Init a Pipeline object.

        Parameters
        ----------
        bids_directory : str, optional
            Path to a BIDS directory. Defaults to None.

        caps_directory : str, optional
            Path to a CAPS directory. Defaults to None.

        tsv_file : str, optional
            Path to a subjects-sessions `.tsv` file. Defaults to None.

        overwrite_caps : bool, optional
            Overwrite or not output directory. Defaults to False.

        base_dir : str, optional
            Working directory (attribute of Nipype::Workflow class). Defaults to None.

        parameters : dict, optional
            Pipeline parameters. Defaults to {}.

        name : str, optional
            Pipeline name. Defaults to None.

        ignore_dependencies : List of str
            List of names of dependencies whose installation checking procedure should be ignored.
            Defaults to None (i.e. all dependencies will be checked).

        Raises
        ------
        RuntimeError: [description]
        """
        import inspect
        import os
        from pathlib import Path
        from tempfile import mkdtemp

        from clinica.utils.caps import (
            build_caps_dataset_description,
            write_caps_dataset_description,
        )
        from clinica.utils.exceptions import ClinicaCAPSError
        from clinica.utils.inputs import check_bids_folder, check_caps_folder

        self._is_built: bool = False
        self._overwrite_caps: bool = overwrite_caps
        self._bids_directory: Optional[Path] = (
            Path(bids_directory).absolute() if bids_directory else None
        )
        self._caps_directory: Optional[Path] = (
            Path(caps_directory).absolute() if caps_directory else None
        )
        self._verbosity: str = "debug"
        self._tsv_file: Optional[Path] = Path(tsv_file).absolute() if tsv_file else None
        self._info_file: Path = (
            Path(os.path.dirname(os.path.abspath(inspect.getfile(self.__class__))))
            / "info.json"
        )
        self._info: Optional[dict] = None
        self._subjects: Optional[List[str]] = None
        self._sessions: Optional[List[str]] = None
        self._input_node = None
        self._output_node = None
        if base_dir:
            self.base_dir = Path(base_dir).absolute()
            self._base_dir_was_specified = True
        else:
            self.base_dir = Path(mkdtemp())
            self._base_dir_was_specified = False

        self._name = name or self.__class__.__name__
        self._parameters = parameters or {}
        self._dependencies: Optional[List[SoftwareDependency]] = None
        self._ignore_dependencies: List[ThirdPartySoftware] = []
        if ignore_dependencies:
            self._ignore_dependencies += [
                ThirdPartySoftware(d) for d in ignore_dependencies
            ]
        self.caps_name = caps_name

        if not self._bids_directory:
            if not self._caps_directory:
                raise RuntimeError(
                    f"The {self._name} pipeline does not contain "
                    "BIDS nor CAPS directory at the initialization."
                )
            try:
                check_caps_folder(self._caps_directory)
            except ClinicaCAPSError as e:
                desc = build_caps_dataset_description(
                    input_dir=self._caps_directory,
                    output_dir=self._caps_directory,
                    processing_name=self._name,
                    dataset_name=self.caps_name,
                    dependencies=self.dependencies,
                )
                raise ClinicaCAPSError(
                    f"{e}\nYou might want to create a 'dataset_description.json' "
                    f"file with the following content:\n{desc}"
                )
            self.is_bids_dir = False
        else:
            check_bids_folder(self._bids_directory)
            self.is_bids_dir = True
            if self._caps_directory is not None:
                if (
                    not self._caps_directory.exists()
                    or len([f for f in self._caps_directory.iterdir()]) == 0
                ):
                    self._caps_directory.mkdir(parents=True, exist_ok=True)
        if self._caps_directory:
            write_caps_dataset_description(
                input_dir=self.input_dir,
                output_dir=self._caps_directory,
                processing_name=self._name,
                dataset_name=self.caps_name,
                dependencies=self.dependencies,
            )
            check_caps_folder(self._caps_directory)
        self._compute_subjects_and_sessions()
        self._init_nodes()

    def _compute_subjects_and_sessions(self):
        from clinica.utils.participant import get_subject_session_list

        self._subjects, self._sessions = get_subject_session_list(
            self.input_dir,
            subject_session_file=self.tsv_file,
            is_bids_dir=self.is_bids_dir,
            use_session_tsv=False,
            tsv_dir=self.base_dir,
        )
        self._subjects, self._sessions = self.filter_qc()

    def filter_qc(self) -> tuple[list[str], list[str]]:
        return self._subjects, self._sessions

    @property
    def input_dir(self) -> Path:
        if self.is_bids_dir:
            return self._bids_directory
        return self._caps_directory

    @property
    def base_dir_was_specified(self) -> bool:
        return self._base_dir_was_specified

    @property
    def is_built(self) -> bool:
        return self._is_built

    @is_built.setter
    def is_built(self, value: bool):
        self._is_built = value

    @property
    def overwrite_caps(self) -> bool:
        return self._overwrite_caps

    @property
    def parameters(self) -> dict:
        return self._parameters

    @parameters.setter
    def parameters(self, value: dict):
        self._parameters = value
        # Need to rebuild input, output and core nodes
        self.is_built = False
        self.init_nodes()

    @property
    def input_node(self) -> Node:
        return self._input_node

    @property
    def output_node(self) -> Node:
        return self._output_node

    @property
    def bids_directory(self) -> Optional[Path]:
        return self._bids_directory

    @property
    def caps_directory(self) -> Optional[Path]:
        return self._caps_directory

    @property
    def subjects(self) -> List[str]:
        return self._subjects

    @subjects.setter
    def subjects(self, value: List[str]):
        self._subjects = value
        self.is_built = False

    @property
    def sessions(self) -> List[str]:
        return self._sessions

    @sessions.setter
    def sessions(self, value: List[str]):
        self._sessions = value
        self.is_built = False

    @property
    def visits(self) -> list[Visit]:
        return [
            Visit(subject, session)
            for subject, session in zip(self.subjects, self.sessions)
        ]

    @visits.setter
    def visits(self, value: list[Visit]):
        self.subjects = [v.subject for v in value]
        self.sessions = [v.session for v in value]

    @property
    def tsv_file(self) -> Optional[Path]:
        return self._tsv_file

    @property
    def info_file(self) -> Path:
        return self._info_file

    def determine_subject_and_session_to_process(self):
        """Query expected output files in the CAPS folder in order to process only those missing.

        If expected output files already exist in the CAPS folder for some subjects and sessions,
        then do not process those again.
        """
        from clinica.utils.stream import log_and_warn

        visits_already_processed = self.get_processed_images()
        if len(visits_already_processed) == 0:
            return
        message = (
            f"In the provided CAPS folder {self.caps_directory}, Clinica found already processed "
            f"images for {len(visits_already_processed)} visit(s):\n- "
        )
        message += "\n- ".join([str(visit) for visit in visits_already_processed])
        message += "\nThose visits will be ignored by Clinica."
        log_and_warn(message, UserWarning)
        self.visits = [
            visit for visit in self.visits if visit not in visits_already_processed
        ]

    @abc.abstractmethod
    def get_processed_images(self) -> list[Visit]:
        """Extract processed image IDs in `caps_directory` based on `subjects`_`sessions`."""
        raise NotImplemented

    def _init_nodes(self) -> None:
        """Init the basic workflow and I/O nodes necessary before build."""
        self._init_input_node()
        self._init_output_node()
        Workflow.__init__(self, self._name, self.base_dir)
        if self.input_node:
            self.add_nodes([self.input_node])
        if self.output_node:
            self.add_nodes([self.output_node])

    def _init_input_node(self) -> None:
        if self.get_input_fields():
            self._input_node = Node(
                name="Input",
                interface=IdentityInterface(
                    fields=self.get_input_fields(), mandatory_inputs=False
                ),
            )

    def _init_output_node(self) -> None:
        if self.get_output_fields():
            self._output_node = Node(
                name="Output",
                interface=IdentityInterface(
                    fields=self.get_output_fields(), mandatory_inputs=False
                ),
            )

    def has_input_connections(self) -> bool:
        """Checks if the Pipeline's input node has been connected.

        Returns
        -------
        bool :
            True if the input node is connected, False otherwise.
        """
        if self.input_node:
            return self._graph.in_degree(self.input_node) > 0
        return False

    def has_output_connections(self) -> bool:
        """Checks if the Pipeline's output node has been connected.

        Returns
        -------
        bool :
            True if the output node is connected, False otherwise.
        """
        if self.output_node:
            return self._graph.out_degree(self.output_node) > 0
        return False

    @postset("is_built", True)
    def build(self):
        """Builds the core, input and output nodes of the Pipeline.

        This method first checks it has already been run. It then checks
        the pipelines dependencies and, in this order, builds the core nodes,
        the input node and, finally, the output node of the Pipeline.

        Since this method returns the concerned object, it can be chained to
        any other method of the Pipeline class.

        Returns
        -------
        self: A Pipeline object.
        """
        if not self.is_built:
            self._check_dependencies()
            self._check_pipeline_parameters()
            if not self.has_input_connections():
                self.determine_subject_and_session_to_process()
                self._build_input_node()
            self._build_core_nodes()
            if not self.has_output_connections():
                self._build_output_node()
        return self

    def run(
        self,
        plugin=None,
        plugin_args=None,
        update_hash: bool = False,
        bypass_check: bool = False,
    ):
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

        from clinica.utils.ux import print_failed_images

        self._handle_cross_sectional_dataset()
        if not self.is_built:
            self.build()
        if not bypass_check:
            self._check_size()
            plugin_args = self._update_parallelize_info(plugin_args)
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
            msg = (
                "Either all the images were already run by the pipeline "
                "or no image was found to run the pipeline."
            )
            log_and_warn(msg, UserWarning)
            exec_graph = Graph()
        return exec_graph

    @property
    def info(self) -> dict:
        import json

        if self._info is None:
            try:
                with open(self.info_file) as info_file:
                    info = json.load(info_file)
            except FileNotFoundError:
                log_and_warn(
                    f"Info file {self.info_file} for pipeline {self.name} is missing.",
                    UserWarning,
                )
                info = {}
            self._info = info
        return self._info

    @property
    def dependencies(self) -> List[SoftwareDependency]:
        if self._dependencies is None:
            dependencies = [
                SoftwareDependency.from_dict(d)
                for d in self.info.get("dependencies", [])
            ]
            dependencies = [
                d for d in dependencies if d.name not in self._ignore_dependencies
            ]
            self._dependencies = dependencies
        return self._dependencies

    def _check_dependencies(self):
        """Checks if listed dependencies are correctly installed."""
        for dependency in self.dependencies:
            dependency.check()
        self._check_custom_dependencies()

    def _check_size(self):
        """Check if the pipeline has enough space on the disk for both working directory and CAPS."""
        import sys
        from os import statvfs
        from os.path import dirname

        from clinica.utils.stream import cprint

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
            space_needed_caps = n_sessions * _human2bytes(space_needed_caps_1_session)
            space_needed_wd = n_sessions * _human2bytes(space_needed_wd_1_session)
            error = ""
            if free_space_caps == free_space_wd:
                if space_needed_caps + space_needed_wd > free_space_wd:
                    # We assume this is the same disk
                    error = (
                        f"{error}Space needed for CAPS and working directory ("
                        f"{_bytes2human(space_needed_caps + space_needed_wd)}"
                        f") is greater than what is left on your hard drive ("
                        f"{_bytes2human(free_space_wd)})"
                    )
            else:
                if space_needed_caps > free_space_caps:
                    error = (
                        f"{error} Space needed for CAPS ("
                        f"{_bytes2human(space_needed_caps)}"
                        f") is greater than what is left on your hard drive ("
                        f"{_bytes2human(free_space_caps)})\n"
                    )
                if space_needed_wd > free_space_wd:
                    error = (
                        f"{error}Space needed for working_directory ("
                        f"{_bytes2human(space_needed_wd)}"
                        f") is greater than what is left on your hard drive ("
                        f"{_bytes2human(free_space_wd)})\n"
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

    def _update_parallelize_info(self, plugin_args: Optional[dict]) -> dict:
        """Performs some checks of the number of threads given in parameters,
        given the number of CPUs of the machine in which clinica is running.
        We force the use of plugin MultiProc

        Author: Arnaud Marcoux
        """
        from multiprocessing import cpu_count

        n_cpu = cpu_count()
        ask_user = False

        try:
            # if no --n_procs arg is used, plugin_arg is None
            # so we need a try / except block
            n_thread_cmdline = plugin_args["n_procs"]
            if n_thread_cmdline > n_cpu:
                msg = (
                    f"You are trying to run clinica with a number of threads ({n_thread_cmdline}) superior to your "
                    f"number of CPUs ({n_cpu})."
                )
                log_and_warn(msg, UserWarning)
                ask_user = True
        except TypeError:
            log_and_warn(
                (
                    f"You did not specify the number of threads to run in parallel (--n_procs argument)."
                    f"Computation time can be shorten as you have {n_cpu} CPUs on this computer. "
                    f"We recommend using {n_cpu - 1} threads."
                ),
                UserWarning,
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

    def _handle_cross_sectional_dataset(self):
        """Check if the dataset is longitudinal or cross-sectional.

        If it is cross-sectional, propose to convert it to a longitudinal layout.
        """
        if self.bids_directory is None:
            return
        subjects = [
            f.name
            for f in self.bids_directory.iterdir()
            if (self.bids_directory / f).is_dir() and f.name.startswith("sub-")
        ]
        (
            cross_sectional_subjects,
            longitudinal_subjects,
        ) = _detect_cross_sectional_and_longitudinal_subjects(
            subjects, self.bids_directory
        )
        if cross_sectional_subjects:
            self._convert_to_longitudinal_if_user_agrees(
                cross_sectional_subjects, longitudinal_subjects
            )

    def _convert_to_longitudinal_if_user_agrees(
        self,
        cross_sectional_subjects: List[str],
        longitudinal_subjects: List[str],
    ):
        import sys

        from clinica.utils.stream import cprint

        log_and_warn(
            (
                f"The following subjects of the input dataset {self.bids_directory} seem to "
                "have a cross-sectional layout which is not supported by Clinica:\n"
                + "\n- ".join(cross_sectional_subjects)
            ),
            UserWarning,
        )
        proposed_bids = (
            self.bids_directory.parent / f"{self.bids_directory.name}_clinica_compliant"
        )
        if not click.confirm(
            "Do you want to proceed with the conversion in another folder? "
            "(Your original BIDS folder will not be modified "
            f"and the folder {proposed_bids} will be created.)"
        ):
            click.echo(
                "Clinica does not support cross-sectional layout for BIDS dataset input. "
                "To run the pipeline, please provide a dataset in longitudinal format or accept "
                "the automatic conversion. Clinica will now exit."
            )
            sys.exit()
        cprint("Converting cross-sectional dataset into longitudinal...")
        _convert_cross_sectional(
            self.bids_directory,
            proposed_bids,
            cross_sectional_subjects,
            longitudinal_subjects,
        )
        cprint(
            f"Conversion succeeded. Your clinica-compliant dataset is located here: {proposed_bids}. "
            "The pipeline will run using this new dataset as input."
        )
        self._bids_directory = proposed_bids
        self._compute_subjects_and_sessions()

    @abc.abstractmethod
    def _build_core_nodes(self):
        """Builds the Pipeline's core nodes.

        This method should use and connect to the `Pipeline.input_node` one
        or more core `Node`s. The outputs of the core processing should then be
        connected to the `Pipeline.output_node`.
        """

    @abc.abstractmethod
    def _build_input_node(self):
        """Builds the Pipeline's input data stream node.

        Warnings:
            This method does not modify the `Pipeline.input_node` (see the
            notes about the global architecture in the class documentation).
        """

    @abc.abstractmethod
    def _build_output_node(self):
        """Builds the Pipeline's output data stream node.

        Warnings:
            This method does not modify the `Pipeline.output_node` (see the
            notes about the global architecture in the class documentation).
        """

    @abc.abstractmethod
    def get_input_fields(self) -> List[str]:
        """Lists the input fields of the Pipeline.

        Returns:
            A list of strings defining the fields of the `IdentityInterface`
            of the `Pipeline.input_node`.
        """

    @abc.abstractmethod
    def get_output_fields(self) -> List[str]:
        """Lists the output fields of the Pipeline.

        Returns:
            A list of strings defining the fields of the `IdentityInterface`
            of the `Pipeline.output_node`.
        """

    @abc.abstractmethod
    def _check_custom_dependencies(self) -> None:
        """Check dependencies provided by the developer."""

    @abc.abstractmethod
    def _check_pipeline_parameters(self) -> None:
        """Check pipeline parameters."""


class GroupPipeline(Pipeline):
    def __init__(
        self,
        caps_directory: str,
        group_label: str,
        bids_directory: Optional[str] = None,
        tsv_file: Optional[str] = None,
        overwrite_caps: Optional[bool] = False,
        base_dir: Optional[str] = None,
        parameters: Optional[dict] = None,
        name: Optional[str] = None,
        ignore_dependencies: Optional[List[str]] = None,
        caps_name: Optional[str] = None,
    ):
        super().__init__(
            bids_directory=bids_directory,
            caps_directory=caps_directory,
            tsv_file=tsv_file,
            overwrite_caps=overwrite_caps,
            base_dir=base_dir,
            parameters=parameters,
            name=name,
            ignore_dependencies=ignore_dependencies,
            caps_name=caps_name,
        )
        self.group_label: GroupLabel = GroupLabel(group_label)
        self.group_id: GroupID = GroupID.from_label(self.group_label)
        self.group_directory: Path = self.caps_directory / "groups" / str(self.group_id)
