"""

"""

import nipype.pipeline.engine as npe
import abc
from bids.grabbids import BIDSLayout


def postset(attribute, value):
    def postset_decorator(func):
        def func_wrapper(self, *args, **kwargs):
            res = func(self, *args, **kwargs)
            setattr(self, attribute, value)
            return res
        return func_wrapper
    return postset_decorator

def get_subject_session_list(input_dir, ss_file=None):

    import clinica.bids.utils.data_handling as cdh
    import pandas as pd

    if not ss_file:
        cdh.create_subs_sess_list(input_dir, '/tmp/')
        ss_file = '/tmp/subjects_sessions_list.tsv'

    ss_df = pd.io.parsers.read_csv(ss_file, sep='\t')
    if list(ss_df.columns.values) != ['participant_id', 'session_id']:
        raise Exception(
            'Subjects and visits file is not in the correct format.')
    subjects = list(ss_df.participant_id)
    sessions = list(ss_df.session_id)

    return sessions, subjects

class Pipeline(npe.Workflow):
    """
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, bids_dir, caps_dir, tsv_file=None):
        """
        """
        self._is_built = False
        self._has_input_node = False
        self._has_output_node = False
        self._input_node = None
        self._output_node = None
        self._input_dir = bids_dir
        self._output_dir = caps_dir
        self._verbosity = 'debug'
        self._tsv_file = tsv_file
        self._name = self.__class__.__name__
        self._parameters = {}
        self._sessions, self._subjects = get_subject_session_list(bids_dir, tsv_file)
        npe.Workflow.__init__(self, self.__class__.__name__)

    def run(self, plugin=None, plugin_args=None, update_hash=False):
        if not self.is_built:
            self.build()
        npe.Workflow.run(self, plugin, plugin_args, update_hash)

    @property
    def is_built(self): return self._is_built

    @is_built.setter
    def is_built(self, value): self._is_built = value

    @property
    def has_input_node(self): return self._has_input_node

    @has_input_node.setter
    def has_input_node(self, value): self._has_input_node = value

    @property
    def has_output_node(self): return self._has_output_node

    @has_output_node.setter
    def has_output_node(self, value): self._has_output_node = value

    @property
    def parameters(self): return self._parameters

    @parameters.setter
    def parameters(self, value): self._parameters = value

    @property
    def input_node(self): return self._input_node

    @input_node.setter
    @postset('has_input_node', True)
    def input_node(self, value): self._input_node = value

    @property
    def output_node(self): return self._output_node

    @output_node.setter
    @postset('has_output_node', True)
    def output_node(self, value): self._output_node = value

    @property
    def input_dir(self): return self._input_dir

    @property
    def output_dir(self): return self._output_dir

    @property
    def subjects(self): return self._subjects

    @property
    def sessions(self): return self._sessions

    @property
    def tsv_file(self): return self._tsv_file

    @property
    def bids_layout(self): return BIDSLayout(self.input_dir)

    @abc.abstractmethod
    def build(self): pass

    @abc.abstractmethod
    def default_input_node(self): pass

    @abc.abstractmethod
    def default_output_node(self): pass

    def check_dependencies(self, dependencies):
        """
        Raise exception if a program in the list do not exist
        Example:
            check_dependencies(["recon-all", "freesurfer"])

        Args:
            dependencies: list of program names

        Raises:
            Exception: Raises an exception.
        """
        from distutils.spawn import find_executable

        def check_dependence(program_name):
            path = find_executable(program_name)
            if path is None:
                raise Exception('Program [%s] do not found' % program_name)

        [check_dependence(x) for x in dependencies]

    @abc.abstractmethod
    def expected_result_files(self): pass
