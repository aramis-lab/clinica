# coding: utf8

# Use hash instead of parameters for iterables folder names
# Otherwise path will be too long and generate OSError
from nipype import config

import clinica.pipelines.engine as cpe

cfg = dict(execution={'parameterize_dirs': False})
config.update_config(cfg)


class T1FreeSurfer(cpe.Pipeline):
    """FreeSurfer-based processing of T1-weighted MR images.

    Attributes:
        subjects: List of participant IDs (e.g. ['sub-CLNC01', 'sub-CLNC01', 'sub-CLNC02', 'sub-CLNC02'])
        sessions: List of session IDs (e.g. ['ses-M00', 'ses-M18', 'ses-M00', 'ses-M18'])
        bids_directory:
        caps_directory:
        recon_all_args: Arguments for recon-all command (e.g. '-qcache')

    Returns:
        A clinica pipeline object containing the T1FreeSurfer pipeline.
    """

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file."""
        pass

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Note:
            The list of inputs of the T1FreeSurfer pipeline is:
                * t1w (str): Path of the T1 weighted image in BIDS format

        Returns:
            A list of (string) input fields name.
        """
        return ['t1w']

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Note:
            The list of outputs of the T1FreeSurfer pipeline is:
                * image_id (str): Image ID (e.g. sub-CLNC01_ses-M00)

        Returns:
            A list of (string) output fields name.
        """
        return ['image_id']

    def build_input_node(self):
        """Build and connect an input node to the pipeline.

        Raise:
            ClinicaBIDSError: If there are duplicated files or missing files for any subject
        """
        import os
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from clinica.utils.exceptions import ClinicaBIDSError, ClinicaException
        from clinica.utils.stream import cprint
        from clinica.utils.ux import (print_images_to_process, print_no_image_to_process)
        from clinica.utils.inputs import clinica_file_reader
        from clinica.utils.input_files import T1W_NII
        from clinica.iotools.utils.data_handling import check_volume_location_in_world_coordinate_system
        from clinica.utils.io import save_participants_sessions

        # Inputs from anat/ folder
        # ========================
        # T1w file:
        try:
            t1w_files = clinica_file_reader(self.subjects,
                                            self.sessions,
                                            self.bids_directory,
                                            T1W_NII)
        except ClinicaException as e:
            err_msg = 'Clinica faced error(s) while trying to read files in your BIDS directory.\n' + str(e)
            raise ClinicaBIDSError(err_msg)

        # Save subjects to process in <WD>/T1FreeSurfer/participants.tsv
        save_participants_sessions(self.subjects, self.sessions, os.path.join(self.base_dir, self.__class__.__name__))

        print_images_to_process(self.subjects, self.sessions)
        cprint('The pipeline will last approximately 10 hours per image.')

        read_node = npe.Node(name="ReadingFiles",
                             iterables=[
                                 ('t1w', t1w_files),
                             ],
                             synchronize=True,
                             interface=nutil.IdentityInterface(
                                 fields=self.get_input_fields())
                             )
        check_volume_location_in_world_coordinate_system(t1w_files, self.bids_directory)
        self.connect([
            (read_node, self.input_node, [('t1w', 't1w')]),
        ])

    def build_output_node(self):
        """Build and connect an output node to the pipeline."""
        import nipype.pipeline.engine as npe
        import nipype.interfaces.utility as nutil
        from .t1_freesurfer_utils import save_to_caps
        import os

        save_to_caps = npe.Node(interface=nutil.Function(
            input_names=['source_dir', 'image_id', 'caps_dir', 'overwrite_caps'],
            output_names=['image_id'],
            function=save_to_caps),
            name='SaveToCaps')
        save_to_caps.inputs.source_dir = os.path.join(self.base_dir, 'T1FreeSurfer', 'ReconAll')
        save_to_caps.inputs.caps_dir = self.caps_directory
        save_to_caps.inputs.overwrite_caps = False

        self.connect([
            (self.output_node, save_to_caps, [('image_id', 'image_id')]),
        ])

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline."""
        import os
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from nipype.interfaces.freesurfer.preprocess import ReconAll
        from .t1_freesurfer_utils import (init_input_node, write_tsv_files)

        # Nodes declaration
        # =================
        # Initialize the pipeline
        #   - Extract <image_id> (e.g. sub-CLNC01_ses-M00) T1w filename;
        #   - Check FOV of T1w;
        #   - Create <subjects_dir> folder in <WD>/T1FreeSurfer/ReconAll/<image_id>/;
        #   - Print begin execution message.
        init_input = npe.Node(
            interface=nutil.Function(
                input_names=['t1w', 'recon_all_args', 'output_dir'],
                output_names=['image_id', 't1w', 'flags', 'subjects_dir'],
                function=init_input_node),
            name='0-InitPipeline')
        init_input.inputs.recon_all_args = self.parameters['recon_all_args']
        init_input.inputs.output_dir = os.path.join(self.base_dir, 'T1FreeSurfer', 'ReconAll')

        # Run recon-all command
        # FreeSurfer segmentation will be in <subjects_dir>/<image_id>/
        recon_all = npe.Node(interface=ReconAll(),
                             name='1-SegmentationReconAll')
        recon_all.inputs.directive = 'all'

        # Generate TSV files containing a summary of the regional statistics
        # in <subjects_dir>/regional_measures
        create_tsv = npe.Node(interface=nutil.Function(
            input_names=['subjects_dir', 'image_id'],
            output_names=['image_id'],
            function=write_tsv_files),
            name='2-CreateTsvFiles')

        # Connections
        # ===========
        self.connect([
            # Get <image_id> from input_node and print begin message
            (self.input_node, init_input, [('t1w', 't1w')]),
            # Run recon-all command
            (init_input, recon_all, [('subjects_dir', 'subjects_dir'),
                                     ('t1w', 'T1_files'),
                                     ('image_id', 'subject_id'),
                                     ('flags',  'flags')]),
            # Generate TSV files
            (init_input, create_tsv, [('subjects_dir', 'subjects_dir')]),
            (recon_all,  create_tsv, [('subject_id', 'image_id')]),
            # Output node
            (recon_all, self.output_node, [('subject_id', 'image_id')]),
        ])
