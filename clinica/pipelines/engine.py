"""

"""

from nipype.pipeline.engine import Workflow
import abc


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


def get_subject_session_list(input_dir, ss_file=None):
    """Parses a BIDS directory to get the subjects and sessions.

    This function lists all the subjects and sessions based on the content of
    the BIDS directory or (if specified) on the provided subject-sessions tsv
    file.

    Args:
        input_dir: A BIDS directory path.
        ss_file: A subjects-sessions file (.tsv format).

    Returns:
        subjects: A subjects list.
        sessions: A sessions list.
    """
    import clinica.iotools.utils.data_handling as cdh
    import pandas as pd

    if not ss_file:
        cdh.create_subs_sess_list(input_dir, '/tmp/') # FIXME(@jguillon): make temporary directory (or temporary file)
        ss_file = '/tmp/subjects_sessions_list.tsv'

    ss_df = pd.io.parsers.read_csv(ss_file, sep='\t')
    if list(ss_df.columns.values) != ['participant_id', 'session_id']:
        raise Exception(
            'Subjects and visits file is not in the correct format.')
    subjects = list(ss_df.participant_id)
    sessions = list(ss_df.session_id)

    return sessions, subjects


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
        bids_layout (:obj:`BIDSLayout`): Object representing the BIDS directory.
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

    def __init__(self, bids_directory=None, caps_directory=None, tsv_file=None, name=None):
        """Inits a Pipeline object.

        Args:
            bids_directory (optional): Path to a BIDS directory.
            caps_directory (optional): Path to a CAPS directory.
            tsv_file (optional): Path to a subjects-sessions `.tsv` file.
            name (optional): A pipelines name.
        """
        import nipype.interfaces.utility as nutil
        import inspect
        import os
        self._is_built = False
        self._bids_directory = bids_directory
        self._caps_directory = caps_directory
        self._verbosity = 'debug'
        self._tsv_file = tsv_file
        self._info_file = os.path.join(
            os.path.dirname(os.path.abspath(inspect.getfile(self.__class__))),
            'info.json')
        self._info = {}
        if name:
            self._name = name
        else:
            self._name = self.__class__.__name__
        self._parameters = {}
        if bids_directory:
            self._sessions, self._subjects = get_subject_session_list(
                bids_directory, tsv_file)
        else:
            self._sessions = []
            self._subjects = []
        if self.get_input_fields():
            self._input_node = npe.Node(name="Input",
                                        interface=nutil.IdentityInterface(
                                            fields=self.get_input_fields(),
                                            mandatory_inputs=False))
        else:
            self._input_node = None
        if self.get_output_fields():
            self._output_node = npe.Node(name="Output",
                                         interface=nutil.IdentityInterface(
                                             fields=self.get_output_fields(),
                                             mandatory_inputs=False))
        else:
            self._output_node = None
        Workflow.__init__(self, self._name)
        if self.input_node: self.add_nodes([self.input_node])
        if self.output_node: self.add_nodes([self.output_node])

    def has_input_connections(self):
        """Checks if the Pipeline's input node has been connected.

        Returns:
            True if the input node is connected, False otherwise.
        """
        return self._graph.in_degree(self.input_node) > 0

    def has_output_connections(self):
        """Checks if the Pipeline's output node has been connected.

        Returns:
            True if the output node is connected, False otherwise.
        """
        return self._graph.out_degree(self.output_node) > 0

    @postset('is_built', True)
    def build(self):
        """Builds the core, input and output nodes of the Pipeline.

        This method first checks it has already been run. It then checks
        the pipelines dependencies and, in this order, builds the core nodes,
        the input node and, finally, the ouput node of the Pipeline.

        Since this method returns the concerned object, it can be chained to
        any other method of the Pipeline class.

        Returns:
            self: A Pipeline object.
        """
        if not self.is_built:
            self.check_dependencies()
            self.build_core_nodes()
            if not self.has_input_connections():
                self.build_input_node()
            if not self.has_output_connections():
                self.build_output_node()
        return self

    def run(self, plugin=None, plugin_args=None, update_hash=False):
        """Executes the Pipeline.

        It overwrites the default Workflow method to check if the
        Pipeline is built before running it. If not, it builds it and then
        run it.

        Args:
            Similar to those of Workflow.run.

        Returns:
            An execution graph (see Workflow.run).
        """
        if not self.is_built:
            self.build()
        return Workflow.run(self, plugin, plugin_args, update_hash)

    def load_info(self):
        """Loads the associated info.json file.

        Todos:
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
    
        Todos:
            - [ ] MATLAB toolbox dependency checking
            - [x] check MATLAB
            - [ ] Clinica pipelines dependency checkings
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
            'ants': chk.check_ants,
            'spm': chk.check_spm,
            'freesurfer': chk.check_freesurfer,
            'fsl': chk.check_fsl,
            'mrtrix': chk.check_mrtrix,
            'matlab': chk.check_matlab
        }
        check_binary = chk.is_binary_present
        # check_toolbox = chk.is_toolbox_present
        # check_pipeline = chk.is_pipeline_present

        # Load the info.json file
        if not self.info:
            self.load_info()

        # Dependencies checking
        for d in self.info['dependencies']:
            if d['type'] == 'software':
                check_software[d['name']]()
            elif d['type'] == 'binary':
                check_binary(d['name'])
            elif d['type'] == 'toolbox':
                pass
            elif d['type'] == 'pipeline':
                pass
            else:
                raise Exception("Unknown dependency type: '%s'." % d['type'])

        self.check_custom_dependencies()

        return self

    @property
    def is_built(self): return self._is_built

    @is_built.setter
    def is_built(self, value): self._is_built = value

    @property
    def parameters(self): return self._parameters

    @parameters.setter
    def parameters(self, value): self._parameters = value

    @property
    def info(self): return self._info

    @info.setter
    def info(self, value): self._info = value

    @property
    def bids_layout(self):
        from bids.grabbids import BIDSLayout
        return BIDSLayout(self.bids_directory)

    @property
    def input_node(self): return self._input_node

    @property
    def output_node(self): return self._output_node

    @property
    def bids_directory(self): return self._bids_directory

    @property
    def caps_directory(self): return self._caps_directory

    @property
    def subjects(self): return self._subjects

    @property
    def sessions(self): return self._sessions

    @property
    def tsv_file(self): return self._tsv_file

    @property
    def info_file(self): return self._info_file

    @abc.abstractmethod
    def build_core_nodes(self):
        """Builds the Pipeline's core nodes.

        This method should use and connect to the `Pipeline.input_node` one
        or more core `Node`s. The outputs of the core processing should then be
        connected to the `Pipeline.output_node`.
        """
        pass

    @abc.abstractmethod
    def build_input_node(self):
        """Builds the Pipeline's input data stream node.

        Warnings:
            This method does not modify the `Pipeline.input_node` (see the
            notes about the global architecture in the class documentation).
        """
        pass

    @abc.abstractmethod
    def build_output_node(self):
        """Builds the Pipeline's output data stream node.

        Warnings:
            This method does not modify the `Pipeline.output_node` (see the
            notes about the global architecture in the class documentation).
        """
        pass

    @abc.abstractmethod
    def get_input_fields(self):
        """Lists the input fields of the Pipeline.

        Returns:
            A list of strings defining the fields of the `IdentityInterface`
            of the `Pipeline.input_node`.
        """
        pass

    @abc.abstractmethod
    def get_output_fields(self):
        """Lists the output fields of the Pipeline.

        Returns:
            A list of strings defining the fields of the `IdentityInterface`
            of the `Pipeline.output_node`.
        """
        pass

    @abc.abstractmethod
    def check_custom_dependencies(self):
        """Checks dependencies provided by the developer.
        """
        pass
