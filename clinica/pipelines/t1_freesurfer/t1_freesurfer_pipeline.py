# coding: utf8

import clinica.pipelines.engine as cpe

# Use hash instead of parameters for iterables folder names
# Otherwise path will be too long and generate OSError
from nipype import config
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
        from colorama import Fore
        from clinica.utils.exceptions import ClinicaBIDSError
        from clinica.utils.io import check_input_bids_files
        from clinica.utils.stream import cprint

        # Remove 'sub-' prefix from participant IDs
        participant_labels = '|'.join(sub[4:] for sub in self.subjects)
        # Remove 'ses-' prefix from session IDs
        session_labels = '|'.join(ses[4:] for ses in self.sessions)

        error_message = ""
        # Inputs from anat/ folder
        # ========================
        # T1w file:
        t1w_files = self.bids_layout.get(type='T1w', return_type='file', extensions=['.nii|.nii.gz'],
                                         subject=participant_labels, session=session_labels)
        error_message += check_input_bids_files(t1w_files, "T1W_NII", self.bids_directory, self.subjects, self.sessions)

        if error_message:
            raise ClinicaBIDSError(error_message)

        if len(t1w_files) == 0:
            import sys
            cprint(
                '%s\nEither all the images were already run by the pipeline or no image was found to run the pipeline. '
                'The program will now exit.%s' % (Fore.BLUE, Fore.RESET))
            sys.exit(0)
        else:
            from clinica.utils.io import save_participants_sessions
            # Save subjects to process in <WD>/T1FreeSurfer/participants.tsv
            save_participants_sessions(self.subjects, self.sessions, os.path.join(self.base_dir, self.__class__.__name__))

            images_to_process = ', '.join(self.subjects[i][4:] + '|' + self.sessions[i][4:]
                                          for i in range(len(self.subjects)))
            cprint('The pipeline will be run on the following image(s): %s' % images_to_process)
            cprint('The pipeline will last approximately 10 hours per image.')

        read_node = npe.Node(name="ReadingFiles",
                             iterables=[
                                 ('t1w', t1w_files),
                             ],
                             synchronize=True,
                             interface=nutil.IdentityInterface(
                                 fields=self.get_input_fields())
                             )
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
        from .t1_freesurfer_utils import (init_input_node, check_flags,
                                          create_subjects_dir, write_tsv_files)
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import os
        from nipype.interfaces.freesurfer.preprocess import ReconAll

        # Nodes declaration
        # =================
        # Get <image_id> (e.g. sub-CLNC01_ses-M00) from input_node
        # and print begin message
        init_input = npe.Node(
            interface=nutil.Function(
                input_names=self.get_input_fields(),
                output_names=['image_id'] + self.get_input_fields(),
                function=init_input_node),
            name='0-InitPipeline')

        # Check recon-all flags given by the user
        # and add -cw256 if needed (i.e. if the FOV of the image is greater than 256)
        check_flags = npe.Node(interface=nutil.Function(
            input_names=['in_t1w', 'recon_all_args'],
            output_names=['flags'],
            function=check_flags),
            name='1-CheckFlag')
        check_flags.inputs.recon_all_args = self.parameters['recon_all_args']

        # Create <subjects_dir> folder for recon-all
        # in <WD>/T1FreeSurfer/ReconAll/<image_id>/
        # Otherwise, recon-all will not run
        create_subjects_dir = npe.Node(interface=nutil.Function(
            input_names=['output_dir', 'image_id'],
            output_names=['subjects_dir'],
            function=create_subjects_dir),
            name='1-CreateSubjectsDir')
        create_subjects_dir.inputs.output_dir = os.path.join(self.base_dir, 'T1FreeSurfer', 'ReconAll')

        # Run recon-all command
        # FreeSurfer segmentation will be in <subjects_dir>/<image_id>/
        # where <subjects_dir> = <WD>/T1FreeSurfer/ReconAll/<image_id>/
        recon_all = npe.Node(interface=ReconAll(),
                             name='2-SegmentationReconAll')
        recon_all.inputs.directive = 'all'

        # Generate TSV files containing a summary of the regional statistics
        # in <subjects_dir>/regional_measures
        # where <subjects_dir> = <WD>/T1FreeSurfer/ReconAll/<image_id>/
        create_tsv = npe.Node(interface=nutil.Function(
            input_names=['subjects_dir', 'image_id'],
            output_names=['image_id'],
            function=write_tsv_files),
            name='3-CreateTsvFiles')

        # Connections
        # ===========
        self.connect([
            # Get <image_id> from input_node and print begin message
            (self.input_node, init_input, [('t1w', 't1w')]),
            # Check recon-all flags given by the user and add -cw256 if needed
            (init_input, check_flags, [('t1w', 'in_t1w')]),
            # Create <subjects_dir> folder for recon-all
            (init_input, create_subjects_dir, [('image_id', 'image_id')]),
            # Run recon-all command
            (check_flags,         recon_all, [('flags',  'flags')]),
            (create_subjects_dir, recon_all, [('subjects_dir', 'subjects_dir')]),
            (init_input,          recon_all, [('image_id', 'subject_id')]),
            (init_input,          recon_all, [('t1w', 'T1_files')]),
            # Generate TSV files
            (create_subjects_dir, create_tsv, [('subjects_dir', 'subjects_dir')]),
            (recon_all,           create_tsv, [('subject_id', 'image_id')]),
            # Output node
            (recon_all, self.output_node, [('subject_id', 'image_id')]),
        ])
