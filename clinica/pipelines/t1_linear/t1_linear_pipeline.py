"""T1 Linear - Clinica Pipeline.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details: https://gitlab.icm-institute.org/aramislab/clinica/wikis/docs/InteractingWithClinica.
"""

# WARNING: Don't put any import statement here except if it's absolutly
# necessary. Put it *inside* the different methods.
# Otherwise it will slow down the dynamic loading of the pipelines list by the
# command line tool.
import clinica.pipelines.engine as cpe


class T1Linear(cpe.Pipeline):
    """T1 Linear SHORT DESCRIPTION.
    This preprocessing pipeline includes globally three steps:
    1) N4 bias correction (performed with ANTS).
    2) Linear registration to MNI (MNI icbm152 nlinear sym template)
      (performed with ANTS) - RegistrationSynQuick.
    3) Cropping the background (in order to save computational power).  4)
    Histogram-based intensity normalization. This is a custom function
    performed by the binary ImageMath included with ANTS.


    Warnings:
        - A WARNING.

    Todos:
        - [x] A FILLED TODO ITEM.
        - [ ] AN ON-GOING TODO ITEM.
    
    Args:
        bids_directory (str): Folder with BIDS structure.  
        caps_directory (str): Folder where CAPS structure will be stored.
        ref_template (str): reference template used for image registration.  
        working_directory (str): Folder containing a temporary space to save
        intermediate results.
        tsv: The Subjects-Sessions list file (in .tsv format).

    Returns:
        A clinica pipeline object containing the T1 Linear pipeline.

    Raises:


    Example:
        >>> from t1_linear import T1Linear
        >>> pipeline = T1Linear('~/MYDATASET_BIDS', '~/MYDATASET_CAPS')
        >>> pipeline.parameters = {
        >>>     # ...
        >>> }
        >>> pipeline.base_dir = '/tmp/'
        >>> pipeline.run()
    """


    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file.
        """
        pass

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Returns:
            A list of (string) input fields name.
        """

        return ['t1w']


    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """

        return ['image_id']


    def build_input_node(self):
        """Build and connect an input node to the pipeline.
        """
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from clinica.utils.inputs import check_bids_folder
        from clinica.utils.participant import get_subject_session_list
        from clinica.utils.exceptions import ClinicaBIDSError, ClinicaException
        from clinica.utils.inputs import clinica_file_reader
        from clinica.utils.input_files import T1W_NII
        from clinica.utils.filemanip import get_subject_id

        check_bids_folder(self.bids_directory)
        is_bids_dir = True

        self.sessions, self.subjects = get_subject_session_list(
                self.bids_directory,
                tsv,
                is_bids_dir,
                False,
                self.base_dir
                )
        # Inputs from anat/ folder
        # ========================
        # T1w file:
        try:
            t1w_files = clinica_file_reader(self.subjects,
                    self.sessions,
                    self.bids_directory,
                    T1W_NII)
        except ClinicaException as e:
            err = 'Clinica faced error(s) while trying to read files in your CAPS directory.\n' + str(e)
            raise ClinicaBIDSError(err)

        # Read tsv file and load inputs
        read_node = npe.Node(name="ReadingFiles",
                iterables=[
                    ('t1w', t1w_files),
                    ],
                synchronize=True,
                interface=nutil.IdentityInterface(
                    fields=get_input_fields())
                )
        
        self.connect([
            (read_node, self.input_node, [('t1w', 't1w')]),
            ])
        
        # This node is supposedly used to load BIDS inputs when this pipeline is
        # not already connected to the output of a previous Clinica pipeline.
        # For the purpose of the example, we simply read input arguments given
        # by the command line interface and transmitted here through the
        # `self.parameters` dictionary and pass it to the `self.input_node` to
        # further by used as input of the core nodes.

    def build_output_node(self):
        """Build and connect an output node to the pipeline.
        """
        # In the same idea as the input node, this output node is supposedly
        # used to write the output fields in a CAPS. It should be executed only
        # if this pipeline output is not already connected to a next Clinica
        # pipeline.
        from nipype.interfaces.io import DataSink
        from nipype.pipeline.ingine as npe


        write_node = npe.Node(
                name="WriteCaps",
                interface=DataSink()
                )
        write_node.inputs.base_directory = caps_directory
        write_node.inputs.parameterization = False
        
        self.connect([
            (self.output_node, write_node, [('image_id', 'image_id')]),
            ])


    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline.
        """

        import t1_linear_utils as utils
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        # Step 1
        # ======
        node1 = npe.Node(name="Step1",
                         interface=nutil.Function(
                             input_names=['in_hello_word'],
                             output_names=[],
                             function=utils.step1))

        # Step 2
        # ======
        node2 = npe.Node(name="Step2",
                         interface=nutil.Function(
                             input_names=['in_hello_word'],
                             output_names=[],
                             function=utils.step2))

        # Connection
        # ==========
        self.connect([
            # STEP 1
            (self.input_node,      node1,    [('hello_word',    'in_hello_word')]),
            # STEP 2
            (self.input_node,      node2,    [('hello_word',    'in_hello_word')])
        ])
