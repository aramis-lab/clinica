# coding: utf8

# WARNING: Don't put any import statement here except if it's absolutely
# necessary. Put it *inside* the different methods.
# Otherwise it will slow down the dynamic loading of the pipelines list by the
# command line tool.
import clinica.pipelines.engine as cpe

# Use hash instead of parameters for iterables folder names
from nipype import config
cfg = dict(execution={'parameterize_dirs': False})
config.update_config(cfg)


class DwiConnectome(cpe.Pipeline):
    """Connectome-based processing of corrected DWI datasets.

    Args:
        input_dir: A BIDS directory.
        output_dir: An empty output directory where CAPS structured data will
        be written.
        subjects_sessions_list: The Subjects-Sessions list file (in .tsv
        format).

    Returns:
        A clinica pipeline object containing the DwiConnectome pipeline.
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

        return ['t1_brain_file', 'wm_mask_file', 'dwi_file',
                'dwi_brainmask_file', 'grad_fsl', 'atlas_files']

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """

        return ['response', 'fod', 'tracts', 'nodes',
                'connectomes']

    def build_input_node(self):
        """Build and connect an input node to the pipeline.
        """
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from clinica.iotools.grabcaps import CAPSLayout
        from clinica.utils.stream import cprint
        from clinica.utils.io import check_input_caps_files
        from clinica.utils.exceptions import ClinicaException, ClinicaCAPSError
        from colorama import Fore
        import re

        # Remove 'sub-' prefix from participant IDs
        participants_regex = '|'.join(sub[4:] for sub in self.subjects)
        # Remove 'ses-' prefix from session IDs
        sessions_regex = '|'.join(ses[4:] for ses in self.sessions)

        # cprint('Loading CAPS folder...')
        caps_layout = CAPSLayout(self.caps_directory)
        # cprint('CAPS folder loaded')

        error_message = ""
        # Inputs from t1-freesurfer pipeline
        # ==================================
        # White matter segmentation
        wm_mask_files = caps_layout.get(freesurfer_file='wm.seg.mgz', return_type='file',
                                        subject=participants_regex, session=sessions_regex)
        error_message += check_input_caps_files(wm_mask_files, "T1_FS_WM", "t1-freesurfer",
                                                self.caps_directory, self.subjects, self.sessions)
        # Desikan parcellation
        aparc_aseg_files = caps_layout.get(freesurfer_file='aparc\+aseg.mgz', return_type='file',
                                           subject=participants_regex, session=sessions_regex)
        error_message += check_input_caps_files(aparc_aseg_files, "T1_FS_DESIKAN", "t1-freesurfer",
                                                self.caps_directory, self.subjects, self.sessions)

        # Destrieux parcellation
        aparc_aseg_a2009s_files = caps_layout.get(freesurfer_file='aparc.a2009s\+aseg.mgz', return_type='file',
                                                  subject=participants_regex, session=sessions_regex)
        error_message += check_input_caps_files(aparc_aseg_a2009s_files, "T1_FS_DESTRIEUX", "t1-freesurfer",
                                                self.caps_directory, self.subjects, self.sessions)

        # Inputs from dwi-preprocessing pipeline
        # ======================================
        # Preprocessed DWI
        dwi_files = caps_layout.get(type='dwi', suffix='preproc', extensions=['.nii|.nii.gz'], return_type='file',
                                    subject=participants_regex, session=sessions_regex)
        error_message += check_input_caps_files(dwi_files, "DWI_PREPROC_NII", "dwi-preprocessing",
                                                self.caps_directory, self.subjects, self.sessions)

        # B0 brainmask
        dwi_brainmask_files = caps_layout.get(type='dwi', suffix='brainmask', extensions=['.nii|.nii.gz'], return_type='file',
                                              subject=participants_regex, session=sessions_regex)
        error_message += check_input_caps_files(dwi_brainmask_files, "DWI_PREPROC_BM", "dwi-preprocessing",
                                                self.caps_directory, self.subjects, self.sessions)

        # Preprocessed bvec
        bvec_files = caps_layout.get(type='dwi', suffix='preproc', return_type='file', extensions='bvec',
                                     subject=participants_regex, session=sessions_regex)
        error_message += check_input_caps_files(bvec_files, "DWI_PREPROC_BVEC", "dwi-preprocessing",
                                                self.caps_directory, self.subjects, self.sessions)

        # Preprocessed bval
        bval_files = caps_layout.get(type='dwi', suffix='preproc', extensions='bval', return_type='file',
                                     subject=participants_regex, session=sessions_regex)
        error_message += check_input_caps_files(bval_files, "DWI_PREPROC_BVAL", "dwi-preprocessing",
                                                self.caps_directory, self.subjects, self.sessions)

        if error_message:
            raise ClinicaCAPSError(error_message)

        # Check space of DWI dataset
        dwi_file_spaces = caps_layout.get(type='dwi', suffix='preproc', extensions=['.nii|.nii.gz'],
                                          target='space', return_type='id',
                                          subject=participants_regex, session=sessions_regex)

        # Return an error if all the DWI files are not in the same space
        if any(a != dwi_file_spaces[0] for a in dwi_file_spaces):
            raise ClinicaCAPSError('Preprocessed DWI files are not all in the '
                                   'same space. Please process them separately '
                                   'using the appropriate subjects/sessions '
                                   '`.tsv` file (-tsv option).')

        # Used only for for T1-B0 registration
        if dwi_file_spaces[0] == 'b0':
            # Brain extracted T1w
            t1_brain_files = caps_layout.get(freesurfer_file='brain.mgz', return_type='file',
                                             subject=participants_regex, session=sessions_regex)
            error_message += check_input_caps_files(t1_brain_files, "T1_FS_BE", "t1-freesurfer",
                                                    self.caps_directory, self.subjects, self.sessions)
            if error_message:
                raise ClinicaCAPSError(error_message)

        list_atlas_files = [
            [aparc_aseg_files[i], aparc_aseg_a2009s_files[i]]
            for i in range(len(self.subjects))
        ]

        list_grad_fsl = [
            (bvec_files[i], bval_files[i])
            for i in range(len(self.subjects))
        ]

        if len(dwi_files) == 0:
            import sys
            cprint('%s\nEither all the images were already run by the pipeline or no image was found to run the pipeline. '
                   'The program will now exit.%s' % (Fore.BLUE, Fore.RESET))
            sys.exit(0)
        else:
            p_id_images_to_process = [re.search(r'(sub-[a-zA-Z0-9]+)', caps_file).group() for caps_file in dwi_files]
            s_id_images_to_process = [re.search(r'(ses-[a-zA-Z0-9]+)', caps_file).group() for caps_file in dwi_files]
            images_to_process = ', '.join(p_id[4:] + '|' + s_id[4:]
                                          for p_id, s_id in zip(p_id_images_to_process, s_id_images_to_process))
            cprint('The pipeline will be run on the following subject(s): %s' % images_to_process)

        if dwi_file_spaces[0] == 'b0':
            self.parameters['dwi_space'] = 'b0'
            read_node = npe.Node(name="ReadingFiles",
                                 iterables=[
                                     ('wm_mask_file', wm_mask_files),
                                     ('t1_brain_file', t1_brain_files),
                                     ('dwi_file', dwi_files),
                                     ('dwi_brainmask_file', dwi_brainmask_files),
                                     ('grad_fsl', list_grad_fsl),
                                     ('atlas_files', list_atlas_files)
                                 ],
                                 synchronize=True,
                                 interface=nutil.IdentityInterface(
                                         fields=self.get_input_fields()))
            self.connect([
                (read_node, self.input_node, [('t1_brain_file', 't1_brain_file')]),
                (read_node, self.input_node, [('wm_mask_file', 'wm_mask_file')]),
                (read_node, self.input_node, [('dwi_file', 'dwi_file')]),
                (read_node, self.input_node, [('dwi_brainmask_file', 'dwi_brainmask_file')]),
                (read_node, self.input_node, [('grad_fsl', 'grad_fsl')]),
                (read_node, self.input_node, [('atlas_files', 'atlas_files')]),
            ])

        elif dwi_file_spaces[0] == 'T1w':
            self.parameters['dwi_space'] = 'T1w'
            read_node = npe.Node(name="ReadingFiles",
                                 iterables=[
                                     ('wm_mask_file', wm_mask_files),
                                     ('dwi_file', dwi_files),
                                     ('dwi_brainmask_file', dwi_brainmask_files),
                                     ('grad_fsl', list_grad_fsl),
                                     ('atlas_files', list_atlas_files)
                                 ],
                                 synchronize=True,
                                 interface=nutil.IdentityInterface(
                                         fields=self.get_input_fields()))
            self.connect([
                (read_node, self.input_node, [('wm_mask_file', 'wm_mask_file')]),
                (read_node, self.input_node, [('dwi_file', 'dwi_file')]),
                (read_node, self.input_node, [('dwi_brainmask_file', 'dwi_brainmask_file')]),
                (read_node, self.input_node, [('grad_fsl', 'grad_fsl')]),
                (read_node, self.input_node, [('atlas_files', 'atlas_files')]),
            ])

        else:
            raise ClinicaCAPSError('Bad preprocessed DWI space. Please check '
                                   'your CAPS folder.')

    def build_output_node(self):
        """Build and connect an output node to the pipeline.
        """
        import nipype.pipeline.engine as npe
        import nipype.interfaces.utility as nutil
        import nipype.interfaces.io as nio
        import clinica.pipelines.dwi_connectome.dwi_connectome_utils as utils

        # Writing CAPS
        # ============
        join_node = npe.JoinNode(name='JoinOutputs',
                                 joinsource='ReadingFiles',
                                 interface=nutil.IdentityInterface(
                                         fields=self.get_output_fields()))

        write_node = npe.MapNode(name='WritingCAPS',
                                 iterfield=['container'] +
                                           ['connectome_based_processing.@' + o
                                            for o in self.get_output_fields()],
                                 interface=nio.DataSink(
                                         infields=['connectome_based_processing.@' + o for o in
                                                   self.get_output_fields()]))
        write_node.inputs.base_directory = self.caps_directory
        write_node.inputs.container = utils.get_containers(self.subjects,
                                                           self.sessions)
        write_node.inputs.substitutions = [('trait_added', '')]
        write_node.inputs.parameterization = False

        self.connect([
            # Writing CAPS
            (self.output_node, join_node, [('response', 'response')]),
            (self.output_node, join_node, [('fod', 'fod')]),
            (self.output_node, join_node, [('tracts', 'tracts')]),
            (self.output_node, join_node, [('nodes', 'nodes')]),
            (self.output_node, join_node, [('connectomes', 'connectomes')]),
            (join_node, write_node, [('response', 'connectome_based_processing.@response')]),
            (join_node, write_node, [('fod', 'connectome_based_processing.@fod')]),
            (join_node, write_node, [('tracts', 'connectome_based_processing.@tracts')]),
            (join_node, write_node, [('nodes', 'connectome_based_processing.@nodes')]),
            (join_node, write_node, [('connectomes', 'connectome_based_processing.@connectomes')]),
        ])

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline.

        Notes:
            - If `FSLOUTPUTTYPE` environment variable is not set, `nipype` takes
            NIFTI by default.

        Todo:
            - [x] Detect space automatically.
            - [ ] Allow for custom parcellations (See TODOs in utils).

        """
        import nipype.interfaces.utility as niu
        import nipype.pipeline.engine as npe
        import nipype.interfaces.fsl as fsl
        import nipype.interfaces.freesurfer as fs
        import nipype.interfaces.mrtrix3 as mrtrix3
        from clinica.lib.nipype.interfaces.mrtrix.preprocess import MRTransform
        from clinica.lib.nipype.interfaces.mrtrix3.reconst import EstimateFOD
        from clinica.lib.nipype.interfaces.mrtrix3.tracking import Tractography
        from clinica.utils.exceptions import ClinicaException, ClinicaCAPSError
        from clinica.utils.stream import cprint
        import clinica.pipelines.dwi_connectome.dwi_connectome_utils as utils
        from clinica.utils.mri_registration import convert_flirt_transformation_to_mrtrix_transformation

        # cprint('Building the pipeline...')

        # Nodes
        # =====

        # B0 Extraction (only if space=b0)
        # -------------
        split_node = npe.Node(name="Reg-0-DWI-B0Extraction",
                              interface=fsl.Split())
        split_node.inputs.output_type = "NIFTI_GZ"
        split_node.inputs.dimension = 't'
        select_node = npe.Node(name="Reg-0-DWI-B0Selection", interface=niu.Select())
        select_node.inputs.index = 0

        # B0 Brain Extraction (only if space=b0)
        # -------------------
        mask_node = npe.Node(name="Reg-0-DWI-BrainMasking",
                             interface=fsl.ApplyMask())
        mask_node.inputs.output_type = "NIFTI_GZ"

        # T1-to-B0 Registration (only if space=b0)
        # ---------------------
        t12b0_reg_node = npe.Node(name="Reg-1-T12B0Registration",
                                  interface=fsl.FLIRT(
                                      dof=6, interp='spline', cost='normmi',
                                      cost_func='normmi',
                                  ))
        t12b0_reg_node.inputs.output_type = "NIFTI_GZ"

        # MGZ File Conversion (only if space=b0)
        # -------------------
        t1_brain_conv_node = npe.Node(name="Reg-0-T1-T1BrainConvertion",
                                      interface=fs.MRIConvert())
        wm_mask_conv_node = npe.Node(name="Reg-0-T1-WMMaskConvertion",
                                     interface=fs.MRIConvert())

        # WM Transformation (only if space=b0)
        # -----------------
        wm_transform_node = npe.Node(name="Reg-2-WMTransformation",
                                     interface=fsl.ApplyXFM())
        wm_transform_node.inputs.apply_xfm = True

        # Nodes Generation
        # ----------------
        label_convert_node = npe.MapNode(name="0-LabelsConversion",
                                         iterfield=['in_file', 'in_config',
                                                    'in_lut', 'out_file'],
                                         interface=mrtrix3.LabelConvert())
        label_convert_node.inputs.in_config = utils.get_conversion_luts()
        label_convert_node.inputs.in_lut = utils.get_luts()

        # FSL flirt matrix to MRtrix matrix Conversion (only if space=b0)
        # --------------------------------------------
        fsl2mrtrix_conv_node = npe.Node(
            name='Reg-2-FSL2MrtrixConversion',
            interface=niu.Function(
                input_names=['in_source_image', 'in_reference_image',
                             'in_flirt_matrix', 'name_output_matrix'],
                output_names=['out_mrtrix_matrix'],
                function=convert_flirt_transformation_to_mrtrix_transformation)
            )

        # Parc. Transformation (only if space=b0)
        # --------------------
        parc_transform_node = npe.MapNode(name="Reg-2-ParcTransformation",
                                          iterfield=["in_files", "out_filename"],
                                          interface=MRTransform())

        # Response Estimation
        # -------------------
        resp_estim_node = npe.Node(name="1a-ResponseEstimation",
                                   interface=mrtrix3.ResponseSD())
        resp_estim_node.inputs.algorithm = 'tournier'

        # FOD Estimation
        # --------------
        fod_estim_node = npe.Node(name="1b-FODEstimation",
                                  interface=EstimateFOD())
        fod_estim_node.inputs.algorithm = 'csd'

        # Tracts Generation
        # -----------------
        tck_gen_node = npe.Node(name="2-TractsGeneration",
                                interface=utils.Tractography())
        tck_gen_node.inputs.n_tracks = self.parameters['n_tracks']
        tck_gen_node.inputs.algorithm = 'iFOD2'

        # BUG: Info package does not exist
        # from nipype.interfaces.mrtrix3.base import Info
        # from distutils.version import LooseVersion
        #
        # if Info.looseversion() >= LooseVersion("3.0"):
        #     tck_gen_node.inputs.select = self.parameters['n_tracks']
        # elif Info.looseversion() <= LooseVersion("0.4"):
        #     tck_gen_node.inputs.n_tracks = self.parameters['n_tracks']
        # else:
        #     from clinica.utils.exceptions import ClinicaException
        #     raise ClinicaException("Your MRtrix version is not supported.")

        # Connectome Generation
        # ---------------------
        # only the parcellation and output filename should be iterable, the tck
        # file stays the same.
        conn_gen_node = npe.MapNode(name="3-ConnectomeGeneration",
                                    iterfield=['in_parc', 'out_file'],
                                    interface=mrtrix3.BuildConnectome())

        # Print begin message
        # -------------------
        print_begin_message = npe.MapNode(
            interface=niu.Function(
                input_names=['in_bids_or_caps_file'],
                function=utils.print_begin_pipeline),
            iterfield='in_bids_or_caps_file',
            name='WriteBeginMessage')

        # Print end message
        # -----------------
        print_end_message = npe.MapNode(
            interface=niu.Function(
                input_names=['in_bids_or_caps_file', 'final_file'],
                function=utils.print_end_pipeline),
            iterfield=['in_bids_or_caps_file'],
            name='WriteEndMessage')

        # CAPS File names Generation
        # --------------------------
        caps_filenames_node = npe.Node(name='CAPSFilenamesGeneration',
                                       interface=niu.Function(
                                           input_names='dwi_file',
                                           output_names=self.get_output_fields(),
                                           function=utils.get_caps_filenames))

        # Connections
        # ===========
        # Computation of the diffusion model, tractography & connectome
        # -------------------------------------------------------------
        self.connect([
            (self.input_node, print_begin_message, [('dwi_file', 'in_bids_or_caps_file')]),  # noqa
            (self.input_node, caps_filenames_node, [('dwi_file', 'dwi_file')]),
            # Response Estimation
            (self.input_node, resp_estim_node, [('dwi_file', 'in_file')]),  # Preproc. DWI # noqa
            (self.input_node, resp_estim_node, [('dwi_brainmask_file', 'in_mask')]),  # B0 brain mask # noqa
            (self.input_node, resp_estim_node, [('grad_fsl', 'grad_fsl')]),  # bvecs and bvals # noqa
            (caps_filenames_node, resp_estim_node, [('response', 'wm_file')]),  # output response filename # noqa
            # FOD Estimation
            (self.input_node, fod_estim_node, [('dwi_file', 'in_file')]),  # Preproc. DWI # noqa
            (resp_estim_node, fod_estim_node, [('wm_file', 'wm_txt')]),  # Response (txt file) # noqa
            (self.input_node, fod_estim_node, [('dwi_brainmask_file', 'mask_file')]),  # B0 brain mask # noqa
            (self.input_node, fod_estim_node, [('grad_fsl', 'grad_fsl')]),  # T1-to-B0 matrix file # noqa
            (caps_filenames_node, fod_estim_node, [('fod', 'wm_odf')]),  # output odf filename # noqa
            # Tracts Generation
            (fod_estim_node, tck_gen_node, [('wm_odf', 'in_file')]),  # ODF file # noqa
            (caps_filenames_node, tck_gen_node, [('tracts', 'out_file')]),  # output tck filename # noqa
            # Label Conversion
            (self.input_node, label_convert_node, [('atlas_files', 'in_file')]),  # atlas image files # noqa
            (caps_filenames_node, label_convert_node, [('nodes', 'out_file')]),  # converted atlas image filenames # noqa
            # Connectomes Generation
            (tck_gen_node,        conn_gen_node, [('out_file', 'in_file')]),  # noqa
            (caps_filenames_node, conn_gen_node, [('connectomes', 'out_file')]),  # noqa
        ])
        # Registration T1-DWI (only if space=b0)
        # -------------------
        if self.parameters['dwi_space'] == 'b0':
            self.connect([
                # MGZ Files Conversion
                (self.input_node, t1_brain_conv_node, [('t1_brain_file', 'in_file')]),  # noqa
                (self.input_node, wm_mask_conv_node, [('wm_mask_file', 'in_file')]),  # noqa
                # B0 Extraction
                (self.input_node, split_node, [('dwi_file', 'in_file')]),  # noqa
                (split_node, select_node, [('out_files', 'inlist')]),  # noqa
                # Masking
                (select_node,     mask_node, [('out', 'in_file')]),  # B0 # noqa
                (self.input_node, mask_node, [('dwi_brainmask_file', 'mask_file')]),  # Brain mask # noqa
                # T1-to-B0 Registration
                (t1_brain_conv_node, t12b0_reg_node, [('out_file', 'in_file')]),  # Brain # noqa
                (mask_node,          t12b0_reg_node, [('out_file', 'reference')]),  # B0 brain-masked # noqa
                # WM Transformation
                (wm_mask_conv_node, wm_transform_node, [('out_file', 'in_file')]),  # Brain mask # noqa
                (mask_node,         wm_transform_node, [('out_file', 'reference')]),  # BO brain-masked # noqa
                (t12b0_reg_node,    wm_transform_node, [('out_matrix_file', 'in_matrix_file')]),  # T1-to-B0 matrix file # noqa
                # FSL flirt matrix to MRtrix matrix Conversion
                (t1_brain_conv_node, fsl2mrtrix_conv_node, [('out_file', 'in_source_image')]),  # noqa
                (mask_node,          fsl2mrtrix_conv_node, [('out_file', 'in_reference_image')]),  # noqa
                (t12b0_reg_node,     fsl2mrtrix_conv_node, [('out_matrix_file', 'in_flirt_matrix')]),  # noqa
                # Apply registration without resampling on parcellations
                (label_convert_node,   parc_transform_node, [('out_file', 'in_files')]),  # noqa
                (fsl2mrtrix_conv_node, parc_transform_node, [('out_mrtrix_matrix', 'linear_transform')]),  # noqa
                (caps_filenames_node,  parc_transform_node, [('nodes', 'out_filename')]),  # noqa
            ])
        # Special care for Parcellation & WM mask
        # ---------------------------------------
        if self.parameters['dwi_space'] == 'b0':
            self.connect([
                (
                wm_transform_node, tck_gen_node, [('out_file', 'seed_image')]),  # noqa
                (parc_transform_node, conn_gen_node, [('out_file', 'in_parc')]),  # noqa
                (parc_transform_node, self.output_node, [('out_file', 'nodes')]),  # noqa
            ])
        elif self.parameters['dwi_space'] == 'T1w':
            self.connect([
                (self.input_node, tck_gen_node, [('wm_mask_file', 'seed_image')]),  # noqa
                (label_convert_node, conn_gen_node, [('out_file', 'in_parc')]),  # noqa
                (label_convert_node, self.output_node, [('out_file', 'nodes')]),  # noqa
            ])
        else:
            raise ClinicaCAPSError(
                    'Bad preprocessed DWI space. Please check your CAPS '
                    'folder.')
        # Outputs
        # -------
        self.connect([
            (resp_estim_node, self.output_node, [('wm_file', 'response')]),
            (fod_estim_node,  self.output_node, [('wm_odf', 'fod')]),
            (tck_gen_node,    self.output_node, [('out_file', 'tracts')]),
            (conn_gen_node,   self.output_node, [('out_file', 'connectomes')]),
            (self.input_node, print_end_message, [('dwi_file', 'in_bids_or_caps_file')]),
            (conn_gen_node,   print_end_message, [('out_file', 'final_file')]),
        ])

        # cprint('Pipeline built')
