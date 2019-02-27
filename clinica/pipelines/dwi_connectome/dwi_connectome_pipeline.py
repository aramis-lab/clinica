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
    """Connectome-based processing of DWI datasets.

    Todos:
        [ ] - Add registration of parcellation onto DWI space

    Args:
        input_dir: A BIDS directory.
        output_dir: An empty output directory where CAPS structured data will
        be written.
        subjects_sessions_list: The Subjects-Sessions list file (in .tsv
        format).

    Returns:
        A clinica pipeline object containing the DwiConnectome pipeline.

    Raises:

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

        from clinica.utils.exceptions import ClinicaCAPSError

        # Reading BIDS files
        # ==================

        # cprint('Loading CAPS folder...')
        caps_layout = CAPSLayout(self.caps_directory)
        # cprint('CAPS folder loaded')

        wm_mask_file = []
        dwi_file_space = []
        t1_brain_file = []
        dwi_file = []
        dwi_brainmask_file = []
        bvec = []
        bval = []
        grad_fsl = []
        aparc_aseg = []
        aparc_aseg_a2009s = []
        atlas_files = []
        cprint('Extracting files...')
        for i in range(len(self.subjects)):
            cprint('\t...subject \'' + str(
                    self.subjects[i][4:]) + '\', session \'' + str(
                    self.sessions[i][4:]) + '\'')

            wm_mask_file.append(caps_layout.get(freesurfer_file='wm.seg.mgz',
                                                session=self.sessions[i][4:],
                                                subject=self.subjects[i][4:],
                                                return_type='file')[0])

            dwi_file_space.append(
                    caps_layout.get(type='dwi', suffix='preproc',
                                    target='space',
                                    session=self.sessions[i][4:],
                                    extensions='nii.gz',
                                    subject=self.subjects[i][4:],
                                    return_type='id')[0])

            if dwi_file_space[i] == 'b0':
                t1_brain_file.append(
                        caps_layout.get(freesurfer_file='brain.mgz',
                                        session=self.sessions[i][4:],
                                        subject=self.subjects[i][4:],
                                        return_type='file')[0])

            aparc_aseg.append(
                    caps_layout.get(freesurfer_file='aparc\+aseg.mgz',
                                    session=self.sessions[i][4:],
                                    subject=self.subjects[i][4:],
                                    return_type='file')[0])

            aparc_aseg_a2009s.append(
                    caps_layout.get(freesurfer_file='aparc.a2009s\+aseg.mgz',
                                    session=self.sessions[i][4:],
                                    subject=self.subjects[i][4:],
                                    return_type='file')[0])

            atlas_files.append([aparc_aseg[i], aparc_aseg_a2009s[i]])

            dwi_file.append(
                    caps_layout.get(type='dwi', suffix='preproc',
                                    space=dwi_file_space[i],
                                    session=self.sessions[i][4:],
                                    extensions='nii.gz',
                                    subject=self.subjects[i][4:],
                                    return_type='file')[0])

            dwi_brainmask_file.append(
                    caps_layout.get(type='dwi', suffix='brainmask',
                                    space=dwi_file_space[i],
                                    session=self.sessions[i][4:],
                                    extensions='nii.gz',
                                    subject=self.subjects[i][4:],
                                    return_type='file')[0])

            bvec.append(
                    caps_layout.get(type='dwi', suffix='preproc',
                                    space=dwi_file_space[i],
                                    session=self.sessions[i][4:],
                                    extensions='bvec',
                                    subject=self.subjects[i][4:],
                                    return_type='file')[0])

            bval.append(
                    caps_layout.get(type='dwi', suffix='preproc',
                                    space=dwi_file_space[i],
                                    session=self.sessions[i][4:],
                                    extensions='bval',
                                    subject=self.subjects[i][4:],
                                    return_type='file')[0])

            grad_fsl.append((bvec[i], bval[i]))

        # Return an error if all the DWI files are not in the same space
        if any(a != dwi_file_space[0] for a in dwi_file_space):
            raise ClinicaCAPSError('Preprocessed DWI files are not all in the '
                                   'same space. Please process them separately '
                                   'using the appropriate subjects/sessions '
                                   '`.tsv` file (-tsv option).')

        elif dwi_file_space[0] == 'b0':
            self.parameters['dwi_space'] = 'b0'
            read_node = npe.Node(name="ReadingFiles",
                                 iterables=[
                                     ('wm_mask_file', wm_mask_file),
                                     ('t1_brain_file', t1_brain_file),
                                     ('dwi_file', dwi_file),
                                     ('dwi_brainmask_file', dwi_brainmask_file),
                                     ('grad_fsl', grad_fsl),
                                     ('atlas_files', atlas_files)
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

        elif dwi_file_space[0] == 'T1w':
            self.parameters['dwi_space'] = 'T1w'
            read_node = npe.Node(name="ReadingFiles",
                                 iterables=[
                                     ('wm_mask_file', wm_mask_file),
                                     ('dwi_file', dwi_file),
                                     ('dwi_brainmask_file', dwi_brainmask_file),
                                     ('grad_fsl', grad_fsl),
                                     ('atlas_files', atlas_files)
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

        Todos:
            - [x] Detect space automatically.
            - [ ] Allow for custom parcellations (See TODOs in utils).

        """
        import nipype.interfaces.utility as niu
        import nipype.pipeline.engine as npe
        import nipype.interfaces.fsl as fsl
        import nipype.interfaces.freesurfer as fs
        import nipype.interfaces.mrtrix3 as mrtrix3
        from clinica.lib.nipype.interfaces.mrtrix3.reconst import EstimateFOD
        from clinica.lib.nipype.interfaces.mrtrix3.tracking import Tractography
        from clinica.utils.exceptions import ClinicaCAPSError
        from clinica.utils.stream import cprint
        import clinica.pipelines.dwi_connectome.dwi_connectome_utils as utils

        # cprint('Building the pipeline...')

        # Nodes
        # =====

        # B0 Extraction
        # -------------
        split_node = npe.Node(name="B0Extraction",
                              interface=fsl.Split())
        split_node.inputs.output_type = "NIFTI_GZ"
        split_node.inputs.dimension = 't'
        select_node = npe.Node(name="B0Selection", interface=niu.Select())
        select_node.inputs.index = 0

        # B0 Brain Extraction
        # -------------------
        mask_node = npe.Node(name="BrainMasking",
                             interface=fsl.ApplyMask())
        mask_node.inputs.output_type = "NIFTI_GZ"

        # T1-to-B0 Registration
        # ---------------------
        t12b0_reg_node = npe.Node(name="T12B0Registration",
                                  interface=fsl.FLIRT())
        t12b0_reg_node.inputs.output_type = "NIFTI_GZ"

        # MGZ File Convertion
        # -------------------
        t1_brain_conv_node = npe.Node(name="T1BrainConvertion",
                                      interface=fs.MRIConvert())
        wm_mask_conv_node = npe.Node(name="WMMaskConvertion",
                                     interface=fs.MRIConvert())

        # WM Transformation
        # -----------------
        wm_transform_node = npe.Node(name="WMTransformation",
                                     interface=fsl.ApplyXFM())
        wm_transform_node.inputs.apply_xfm = True

        # Response Estimation
        # -------------------
        resp_estim_node = npe.Node(name="ResponseEstimation",
                                   interface=mrtrix3.ResponseSD())
        resp_estim_node.inputs.algorithm = 'tournier'

        # FOD Estimation
        # --------------
        fod_estim_node = npe.Node(name="FODEstimation",
                                  interface=EstimateFOD())
        fod_estim_node.inputs.algorithm = 'csd'

        # Tracts Generation
        # -----------------
        tck_gen_node = npe.Node(name="TractsGeneration",
                                interface=utils.Tractography())
        tck_gen_node.inputs.n_tracks = self.parameters['n_tracks']
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

        # Nodes Generation
        # ----------------
        label_convert_node = npe.MapNode(name="LabelsConversion",
                                         iterfield=['in_file', 'in_config',
                                                    'in_lut', 'out_file'],
                                         interface=mrtrix3.LabelConvert())
        label_convert_node.inputs.in_config = utils.get_conversion_luts()
        label_convert_node.inputs.in_lut = utils.get_luts()

        # Connectome Generation
        # ---------------------
        # only the parcellation and output filename should be iterable, the tck
        # file stays the same.
        conn_gen_node = npe.MapNode(name="ConnectomeGeneration",
                                    iterfield=['in_parc', 'out_file'],
                                    interface=mrtrix3.BuildConnectome())

        # Print begin message
        # -------------------
        print_begin_message = npe.Node(
            interface=niu.Function(
                input_names=['in_bids_or_caps_file'],
                function=utils.print_begin_pipeline),
            iterfield='in_bids_or_caps_file',
            name='WriteBeginMessage')

        # Print end message
        # -----------------
        print_end_message = npe.Node(
            interface=niu.Function(
                input_names=['in_bids_or_caps_file', 'final_file', 'final_file'],
                function=utils.print_end_pipeline),
            iterfield=['in_bids_or_caps_file', 'final_file'],
            name='WriteEndMessage')

        # CAPS FIlenames Generation
        # -------------------------
        caps_filenames_node = npe.Node(name='CAPSFilenamesGeneration',
                                       interface=niu.Function(
                                               input_names='dwi_file',
                                               output_names=self.get_output_fields(),
                                               function=utils.get_caps_filenames))

        # Connections
        # ===========
        self.connect([
            (self.input_node, print_begin_message, [('dwi_file', 'in_bids_or_caps_file')]),  # noqa
        ])
        # WM in B0-space Transformation
        # -----------------------------
        if self.parameters['dwi_space'] == 'b0':
            self.connect([
                # MGZ Files Convertion
                (self.input_node, t1_brain_conv_node,
                 [('t1_brain_file', 'in_file')]),

                # noqa
                (self.input_node, wm_mask_conv_node,
                 [('wm_mask_file', 'in_file')]),
                # noqa
                # B0 Extraction
                (self.input_node, split_node, [('dwi_file', 'in_file')]),
                # noqa
                (split_node, select_node, [('out_files', 'inlist')]),
                # noqa
                # Masking
                (select_node, mask_node, [('out', 'in_file')]),  # B0 #noqa
                (self.input_node, mask_node,
                 [('dwi_brainmask_file', 'mask_file')]),
                # Brain mask #noqa #noqa
                # T1-to-B0 Registration
                (t1_brain_conv_node, t12b0_reg_node, [('out_file', 'in_file')]),
                # Brain #noqa
                (mask_node, t12b0_reg_node, [('out_file', 'reference')]),
                # B0 brain-masked #noqa
                # WM Transformation
                (wm_mask_conv_node, wm_transform_node,
                 [('out_file', 'in_file')]),  # Brain mask #noqa
                (mask_node, wm_transform_node, [('out_file', 'reference')]),
                # BO brain-masked #noqa
                (t12b0_reg_node, wm_transform_node,
                 [('out_matrix_file', 'in_matrix_file')]),
                # T1-to-B0 matrix file #noqa
            ])

        # Tractography
        # ------------
        self.connect([
            (self.input_node, caps_filenames_node, [('dwi_file', 'dwi_file')]),
            # Response Estimation
            (self.input_node, resp_estim_node, [('dwi_file', 'in_file')]),
            # Preproc. DWI #noqa
            (self.input_node, resp_estim_node,
             [('dwi_brainmask_file', 'in_mask')]),  # B0 brain mask #noqa
            (self.input_node, resp_estim_node, [('grad_fsl', 'grad_fsl')]),
            # bvecs and bvals #noqa
            (caps_filenames_node, resp_estim_node, [('response', 'wm_file')]),
            # output response filename #noqa
            # FOD Estimation
            (self.input_node, fod_estim_node, [('dwi_file', 'in_file')]),
            # Preproc. DWI #noqa
            (resp_estim_node, fod_estim_node, [('wm_file', 'wm_txt')]),
            # Response (txt file) #noqa
            (self.input_node, fod_estim_node,
             [('dwi_brainmask_file', 'mask_file')]),  # B0 brain mask #noqa
            (self.input_node, fod_estim_node, [('grad_fsl', 'grad_fsl')]),
            # T1-to-B0 matrix file #noqa
            (caps_filenames_node, fod_estim_node, [('fod', 'wm_odf')]),
            # output odf filename #noqa
            # Tracts Generation
            (fod_estim_node, tck_gen_node, [('wm_odf', 'in_file')]),
            # ODF file #noqa
            (caps_filenames_node, tck_gen_node, [('tracts', 'out_file')]),
            # output tck filename #noqa
            # Label Conversion
            (self.input_node, label_convert_node, [('atlas_files', 'in_file')]),
            # atlas image files #noqa
            (caps_filenames_node, label_convert_node, [('nodes', 'out_file')]),
            # converted atlas image filenames #noqa
            # Connectomes Generation
            (tck_gen_node, conn_gen_node, [('out_file', 'in_file')]),
            # output odf filename #noqa
            (label_convert_node, conn_gen_node, [('out_file', 'in_parc')]),
            # output odf filename #noqa
            (caps_filenames_node, conn_gen_node, [('connectomes', 'out_file')]),
            # output odf filename #noqa
        ])

        if self.parameters['dwi_space'] == 'b0':
            self.connect([
                (wm_transform_node, tck_gen_node, [('out_file', 'seed_image')])
                # ODF file #noqa
            ])
        elif self.parameters['dwi_space'] == 'T1w':
            self.connect([
                (
                    self.input_node, tck_gen_node,
                    [('wm_mask_file', 'seed_image')])
                # ODF file #noqa
            ])
        else:
            raise ClinicaCAPSError(
                    'Bad preprocessed DWI space. Please check your CAPS '
                    'folder.')

        # Outputs
        # ------
        self.connect([
            (resp_estim_node, self.output_node, [('wm_file', 'response')]),
            # T1-to-B0 matrix file #noqa
            (fod_estim_node, self.output_node, [('wm_odf', 'fod')]),
            # T1-to-B0 matrix file #noqa
            (tck_gen_node, self.output_node, [('out_file', 'tracts')]),
            # T1-to-B0 matrix file #noqa
            (label_convert_node, self.output_node, [('out_file', 'nodes')]),
            # T1-to-B0 matrix file #noqa
            (conn_gen_node, self.output_node, [('out_file', 'connectomes')]),
            # T1-to-B0 matrix file #noqa
            (self.input_node, print_end_message, [('dwi_file', 'in_bids_or_caps_file')]),  # noqa
            (conn_gen_node,   print_end_message, [('out_file', 'final_file')]),  # noqa
        ])

        # cprint('Pipeline built')
