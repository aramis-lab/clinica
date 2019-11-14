# -*- coding: utf-8 -*-

"""deformetrica_preprocessing - Clinica Pipeline.
"""

import clinica.pipelines.engine as cpe

__author__ = "Alexis Guyot"
__copyright__ = "Copyright 2016-2019, The Aramis Lab Team"
__credits__ = ["Alexis Guyot"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Alexis Guyot"
__email__ = "alexis.guyot@icm-institute.org"
__status__ = "Development"


class DeformetricaPrepareData(cpe.Pipeline):
    """Preprocessing for Deformetrica software

    Creates a pipeline that extract any brain structure from a set of
    T1 MR images, align those in the same reference space and generate a
    template from them using deformetrica.

    Args:
        input_dir (String): A BIDS directory.
        output_dir (String): An empty output directory where CAPS
            structured data will be written.
        subjects_sessions_list: The Subjects-Sessions list file (in .tsv
            format).

    Returns:
        A clinica pipeline object containing the DeformetricaPrepareData
            pipeline.
    """
    def __init__(
            self,
            bids_directory=None,
            caps_directory=None,
            tsv_file=None,
            name=None,
            group_id='default'):
        import os

        super(DeformetricaPrepareData, self).__init__(
            bids_directory,
            caps_directory,
            tsv_file,
            name)

        if not group_id.isalnum():
            raise ValueError(
                'Not a valid group_id value. It must be composed only of letters and/or numbers')

        # Check that group does not already exists
        path_exists = os.path.exists(
            os.path.join(
                os.path.abspath(caps_directory),
                'groups',
                'group-{0}'.format(group_id)))
        if path_exists:
            error_message = 'group_id: {0} already exists, please choose another one. Groups that exist in your CAPS directory are: \n'.format(group_id)
            list_groups = os.listdir(
                os.path.join(
                    os.path.abspath(caps_directory),
                    'groups'))
            for e in list_groups:
                if e.startswith('group-'):
                    error_message += '{0} \n'.format(e)
            raise ValueError(error_message)

        # Check that there are at least 2 subjects
        if len(self.subjects) <= 1:
            raise ValueError(
                str.format(
                    'This pipeline needs at least 2 subjects to produce input to [Deformetrica - generate template], and found {0} only in {1}.',
                    len(self.subjects),
                    self.tsv_file))

        self._group_id = group_id

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in `info.json`.
        """
        import nipype.interfaces.vtkbase as vtkbase

        # vtk
        vtk_version = vtkbase.vtk_version()
        if not vtk_version:
            raise RuntimeError('No version of VTK found on the system.')

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Returns:
            A list of (string) input fields name.
        """
        return [
            'input_brains',
            'input_segmentations']

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """
        return []

    def build_input_node(self):
        """Build and connect an input node to the pipeline.
        """
        import os
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import clinica.pipelines.deformetrica_prepare_data.deformetrica_prepare_data_utils as utils
        from clinica.utils.exceptions import ClinicaBIDSError
        from clinica.utils.inputs import clinica_file_reader
        from clinica.utils.input_files import T1_FS_BRAIN

        # An IdentityInterface to get the T1 CAPS segmentations and
        # brain masks
        try:
            caps_brains = clinica_file_reader(self.subjects,
                                              self.sessions,
                                              self.caps_directory,
                                              T1_FS_BRAIN)
            T1_FS_ASEG = {
                'pattern': 't1/freesurfer_cross_sectional/sub-*_ses-*/mri/aseg.mgz',
                'description': ' extracted segmentations from T1w MRI',
                'needed_pipeline': 't1-freesurfer'}
            caps_segmentations = clinica_file_reader(self.subjects,
                                                     self.sessions,
                                                     self.caps_directory,
                                                     T1_FS_ASEG)
        except ClinicaException as e:
            err_msg = 'Clinica faced error(s) while trying to read files in your BIDS directory.\n' + str(e)
            raise ClinicaBIDSError(err_msg)

        read_caps_node = npe.Node(name="read_segmentations_node",
                                  interface=nutil.IdentityInterface(
                                      fields=self.get_input_fields()))
        read_caps_node.inputs.input_brains = caps_brains
        read_caps_node.inputs.input_segmentations = caps_segmentations


        self.connect([
            (
                read_caps_node, self.input_node,
                [('input_brains', 'input_brains')]),
            (
                read_caps_node, self.input_node,
                [('input_segmentations', 'input_segmentations')])])

    def build_output_node(self):
        """Build and connect an output node to the pipeline.

        Similarly to the input node, this output node is supposedly
        used to write the output fields in a CAPS. It should be executed
        only if this pipeline output is not already connected to a next
        Clinica pipeline.
        """
        pass

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline.
        """
        import os

        import clinica.pipelines.deformetrica_prepare_data.deformetrica_prepare_data_utils as utils
        import nipype.interfaces.io as nio
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.freesurfer as freesurfer
        import nipype.interfaces.fsl as fsl

        # Step 1: Convert brain masks to .nii format (for FSL flirt)
        # ======
        brain_nii_node_name = 'step1_brain_nii_node'
        brain_nii_node = npe.MapNode(
            name=brain_nii_node_name,
            interface=freesurfer.MRIConvert(),
            iterfield=['in_file'])
        brain_nii_node.inputs.out_type = 'nii'
        # force to same orientation as Colin 27 volume
        brain_nii_node.inputs.out_orientation = 'RAS'

#        # Step 2: get affine matrix T1 (brain mask) list <-> Colin 27
#        # ======
#        get_affine_node_name = 'step2_get_affine_node'
#        get_affine_node = npe.MapNode(
#            name=get_affine_node_name,
#            interface=fsl.FLIRT(),
#            iterfield=['in_file'])
#        colin27_image_path = os.path.abspath(
#            os.path.join(
#                os.path.dirname(os.path.realpath(__file__)),
#                '..',
#                '..',
#                'resources',
#                'mni',
#                'colin27',
#                'colin27_t1_tal_lin_brain.nii.gz'))
#        get_affine_node.inputs.reference = utils.select_colin27_image(
#            colin27_image_path)
#
#        # Step 3: process [structure]/[structure file] input
#        # ======
#        process_structure_node_name = 'step3_process_structure_node'
#        process_structure_node = npe.Node(
#            name=process_structure_node_name,
#            interface=nutil.Function(
#                input_names=[
#                    'in_structure',
#                    'in_structure_file',
#                    'in_colin27_resources_folder'],
#                output_names=[
#                    'out_structure_list',
#                    'out_structure_file'],
#                function=utils.process_structure))
#        process_structure_node.inputs.in_structure = self.parameters['structure']
#        process_structure_node.inputs.in_structure_file = self.parameters['structure_file']
#        colin27_resources_path = os.path.abspath(
#            os.path.join(
#                os.path.dirname(os.path.realpath(__file__)),
#                '..',
#                '..',
#                'resources',
#                'mni',
#                'colin27'))
#        process_structure_node.inputs.in_colin27_resources_folder = colin27_resources_path
#
#        # Step 4: Get Colin 27 structure templates
#        # =======
#        get_colin27_structure_template_node_name = 'step4_get_colin27_structure_path_node'
#        get_colin27_structure_template_node = npe.MapNode(
#            name=get_colin27_structure_template_node_name,
#            interface=nutil.Function(
#                input_names=[
#                    'in_structure',
#                    'in_colin27_resources_folder'],
#                output_names=[
#                    'out_colin27_structure_template'],
#                function=utils.get_colin27_structure_template
#                ),
#            iterfield=['in_structure'])
#        get_colin27_structure_template_node.inputs.in_colin27_resources_folder = colin27_resources_path
#
#        # Step 5: Link several objects to structures so we can loop over
#        #         ([object], structure)
#        # ======
#        # input: list of objects, list of structures
#        # output: list of all combinations {object, structure} + map
#        #         linking each object to all corresponding
#        #         {object, structure}
#        link_objects_structures_node_name = 'step5_link_objects_structures_node'
#        link_objects_structures_node = npe.Node(
#            name=link_objects_structures_node_name,
#            interface=nutil.Function(
#                input_names=[
#                    'in_brain_list',
#                    'in_segmentation_list',
#                    'in_affine_list',
#                    'in_structure_list'],
#                output_names=[
#                    'out_brain_list',
#                    'out_segmentation_list',
#                    'out_affine_list',
#                    'out_structure_list',
#                    'out_map'],
#                function=utils.link_objects_structures))
#
#        # Step 6: Get list of structure IDs from list of structures
#        # ======
#        get_structureids_node_name = 'step6_get_structureids_node'
#        get_structureids_node = npe.MapNode(
#            name=get_structureids_node_name,
#            interface=nutil.Function(
#                input_names=[
#                    'in_structure'],
#                output_names=[
#                    'out_structure_id'],
#                function=utils.get_structure_id
#                ),
#            iterfield=['in_structure'])
#
#        # Step 7: Tesselate structure from T1 recon-all output
#        # ======
#        # Note: we're iterating over each {subject, structure}
#        tesselate_node_name = 'step7_tesselate_node'
#        tesselate_node = npe.MapNode(
#            name=tesselate_node_name,
#            interface=freesurfer.MRITessellate(),
#            iterfield=['in_file', 'label_value'])
#        # Force to RAS coordinates so both subject surface and Colin27
#        # surface are in RAS
#        tesselate_node.inputs.use_real_RAS_coordinates = True
#
#        # Step 8: Smooth structure
#        # ======
#        smooth_node_name = 'step8_smooth_node'
#        smooth_node = npe.MapNode(
#            name=smooth_node_name,
#            interface=freesurfer.SmoothTessellation(),
#            iterfield=['in_file'])
#
#        # Step 9: Return .vtk structure surfaces
#        # ======
#        get_vtk_surface_node_name = 'step9_get_vtk_surface_node'
#        get_vtk_surface_node = npe.MapNode(
#            name=get_vtk_surface_node_name,
#            interface=freesurfer.MRIsConvert(),
#            iterfield=['in_file'])
#        get_vtk_surface_node.inputs.out_datatype = 'vtk'
#
#        # Step 10: Apply affine transformation to meshes
#        # =======
#        apply_affine_node_name = 'step10_apply_affine_node'
#        apply_affine_node = npe.MapNode(
#            name=apply_affine_node_name,
#            interface=nutil.Function(
#                input_names=[
#                    'in_segmentation',
#                    'in_affine',
#                    'in_brain',
#                    'in_colin27',
#                    'in_structure'],
#                output_names=[
#                    'out_affine_reoriented'],
#                function=utils.apply_affine),
#            iterfield=[
#                'in_segmentation',
#                'in_affine',
#                'in_brain',
#                'in_structure'])
#        get_affine_node.inputs.in_colin27 = utils.select_colin27_image(
#            colin27_image_path)
#
#        # Step 11: get per-subject list of meshes
#        # =======
#        # Note: For each subject, we want the list of mesh
#        #     corresponding to
#        #     structure 1, structure 2, ..., structure n
#        # takes as input a map from {subject} list + {structure} list
#        # to {subject, structure} list
#        get_persubject_meshes_node_name = 'step11_get_persubject_meshes_node'
#        get_persubject_meshes_node = npe.Node(
#            name=get_persubject_meshes_node_name,
#            interface=nutil.Function(
#                input_names=[
#                    'in_mesh_list',
#                    'in_map'],
#                output_names=[
#                    'out_persubject_mesh_list'],
#                function=utils.get_persubject_meshes))
#
#        # Step 12: Rigid reg subject structure mesh <-> Colin 27
#        #          structure mesh
#        # =======
#        rigid_reg_node_name = 'step12_rigid_reg_node'
#        rigid_reg_node = npe.MapNode(
#            name=rigid_reg_node_name,
#            interface=nutil.Function(
#                input_names=[
#                    'in_persubject_mesh_list',
#                    'in_colin27_mesh',
#                    'in_structure_list'],
#                output_names=[
#                    'out_persubject_regmesh_list'],
#                function=utils.rigid_register_colin27
#                ),
#            iterfield=['in_persubject_mesh_list'])
#
#        # Step 13: create per-{subject-structure} datasink basedir
#        # =======
#        # This links the persubject datasink base directory to the list
#        # of structures
#        link_basedir_structures_node_name = 'step13_link_basedir_structures_node'
#        link_basedir_structures_node = npe.Node(
#            name=link_basedir_structures_node_name,
#            interface=nutil.Function(
#                input_names=[
#                    'in_basedir_list',
#                    'in_structure_list'],
#                output_names=[
#                    'out_basedir_list'],
#                function=utils.link_basedirs_structures))
#        subject_data_basedir_list = [
#            os.path.join(
#                os.path.abspath(self.caps_directory),
#                'subjects',
#                subj_sess[0], # subject
#                subj_sess[1], # session
#                't1') for subj_sess in zip(self.subjects, self.sessions)]
#        link_basedir_structures_node.inputs.in_basedir_list = subject_data_basedir_list
#
#        # Step 14: save per-{subject-structure} data to CAPS_DIR
#        # =======
#        # - initial .vtk meshes
#        # - affine registered .vtk meshes
#        subject_structure_datasink_name = 'step14_subject_structure_datasink'
#        data_subfolder = "subfolder"
#        subject_structure_datasink = npe.MapNode(
#            name=subject_structure_datasink_name,
#            interface=nio.DataSink(
#                infields=[
#                    '{0}.@initial_surface'.format(data_subfolder),
#                    '{0}.@affine_surface'.format(data_subfolder)]),
#            iterfield=[
#                'base_directory',
#                '{0}.@initial_surface'.format(data_subfolder),
#                '{0}.@affine_surface'.format(data_subfolder)])
#        subject_structure_data_regexp_substitution_list = utils.initialise_mapnode_datasink_regexp(
#            subject_structure_datasink_name)
#        # Rename affine registered surfaces
#        affine_surface_substitution = (
#            str.format(
#                r'sub-(.*)/ses-(.*)/t1/{0}/_{1}.*/affine_registered_(.*)\.vtk',
#                data_subfolder,
#                apply_affine_node_name),
#            r'sub-\1/ses-\2/t1/target-MNIColin27_AffineRegistration/sub-\1_ses-\2_target-MNIColin27_affine_\3.vtk')
#        subject_structure_data_regexp_substitution_list.append(
#            affine_surface_substitution)
#        # Rename initial surfaces for all potential structures
#        valid_structure_list = utils.get_valid_structure_list(colin27_resources_path)
#        for valid_structure in valid_structure_list:
#            valid_structure_id = utils.get_structure_id(valid_structure)
#            # rename initial surface given the structure name / corresponding
#            # Freesurfer ID
#            initial_substitution = (
#                str.format(
#                    r'sub-(.*)/ses-(.*)/t1/{0}/_{1}.*/aseg_smoothed\.mgz_{2}_converted\.vtk',
#                    data_subfolder,
#                    get_vtk_surface_node,
#                    valid_structure_id),
#                str.format(
#                    r'sub-\1/ses-\2/t1/freesurfer_cross_sectional/sub-\1_ses-\2/bem/sub-\1_ses-\2_orientation-RAS_{0}.vtk',
#                    valid_structure))
#            subject_structure_data_regexp_substitution_list.append(
#                initial_substitution)
#        subject_structure_datasink.inputs.regexp_substitutions = subject_structure_data_regexp_substitution_list
#
#        # Step 15: save per-{subject} data to CAPS_DIR
#        # =======
#        # - rigidly registered .vtk meshes
#        subject_datasink_name = 'step15_subject_datasink'
#        subject_datasink = npe.MapNode(
#            name=subject_datasink_name,
#            interface=nio.DataSink(
#                infields=[
#                    '{0}.@rigid_surface'.format(data_subfolder),
#                    '{0}.@affine_matrix'.format(data_subfolder),
#                    '{0}.@affine_volume'.format(data_subfolder)]),
#            iterfield=[
#                'base_directory',
#                '{0}.@rigid_surface'.format(data_subfolder),
#                '{0}.@affine_matrix'.format(data_subfolder),
#                '{0}.@affine_volume'.format(data_subfolder)])
#        subject_datasink.inputs.base_directory = subject_data_basedir_list
#        subject_data_regexp_substitution_list = utils.initialise_mapnode_datasink_regexp(
#            subject_datasink_name)
#        # rename rigidly registered surfaces
#        rigid_surface_substitution = (
#            str.format(
#                r'sub-(.*)/ses-(.*)/t1/{0}/_{1}.*/out_reoriented_(.*).vtk',
#                data_subfolder,
#                rigid_reg_node_name),
#            r'sub-\1/ses-\2/t1/target-MNIColin27_RigidRegistration/sub-\1_ses-\2_target-MNIColin27_rigid_\3.vtk') 
#        # rename affine registration matrices
#        affine_matrix_substitution = (
#            str.format(
#                r'sub-(.*)/ses-(.*)/t1/{0}/_{1}.*/brainmask_out_flirt.mat',
#                data_subfolder,
#                get_affine_node_name),
#            r'sub-\1/ses-\2/t1/target-MNIColin27_AffineRegistration/sub-\1_ses-\2_target-MNIColin27_affine-fsl.mat') 
#        # rename affine registration matrices
#        affine_volume_substitution = (
#            str.format(
#                r'sub-(.*)/ses-(.*)/t1/{0}/_{1}.*/brainmask_out_flirt.nii.gz',
#                data_subfolder,
#                get_affine_node_name),
#            r'sub-\1/ses-\2/t1/target-MNIColin27_AffineRegistration/sub-\1_ses-\2_target-MNIColin27_T1w.nii.gz') 
#        subject_data_regexp_substitution_list.append(
#            rigid_surface_substitution)
#        subject_data_regexp_substitution_list.append(
#            affine_matrix_substitution)
#        subject_data_regexp_substitution_list.append(
#            affine_volume_substitution)
#        subject_datasink.inputs.regexp_substitutions = subject_data_regexp_substitution_list
#
#        # Step 16: save per-{structure} data to CAPS_DIR
#        # =======
#        # - Colin 16 structure templates (paths)
#        structure_datasink_name = 'step16_structure_datasink'
#        structure_datasink = npe.MapNode(
#            name=structure_datasink_name,
#            interface=nio.DataSink(infields=['templates']),
#            iterfield=['templates'])
#        structure_datasink.inputs.base_directory = os.path.join(
#            self.caps_directory,
#            'groups',
#            'group-{0}'.format(self._group_id))
#        structure_data_regexp_substitution_list = utils.initialise_mapnode_datasink_regexp(
#            structure_datasink_name)
#        # add 'template' to structure name
#        structure_data_regexp_substitution_list.append((
#            str.format(
#                r'groups/group-{0}/templates/colin27_(.*).vtk',
#                self._group_id),
#            str.format(
#                r'groups/group-{0}/templates/colin27_template_\1.vtk',
#                self._group_id)))
#        structure_datasink.inputs.regexp_substitutions = structure_data_regexp_substitution_list
#
#        # Step 17: save common group data to CAPS_DIR
#        # =======
#        # - Deformetrica parameter files (optimisation, data, model)
#        # - initial .vtk template
#        # - Colin 27 .nii.gz brain volume
#        group_datasink_name = 'step17_group_datasink'
#        group_datasink = npe.Node(
#            name=group_datasink_name,
#            interface=nio.DataSink())
#        group_datasink.inputs.base_directory = os.path.join(
#            self.caps_directory,
#            'groups',
#            'group-{0}'.format(self._group_id))
#        group_data_regexp_substitution_list = [(
#            str.format(
#                r'groups/group-{0}/structures/out_structure_file.txt',
#                self._group_id),
#            str.format(
#                r'groups/group-{0}/structures/structure_list.txt',
#                self._group_id))]
#        group_datasink.inputs.regexp_substitutions = group_data_regexp_substitution_list

        # Connection
        # ==========
        self.connect([
            # .nii conversion of brain masks
            (
                self.input_node, brain_nii_node,
                [('input_brains', 'in_file')]),
#            # affine registration
#            (
#                brain_nii_node, get_affine_node,
#                [('out_file', 'in_file')]),
#            # get path to Colin 27 structures
#            (
#                process_structure_node, get_colin27_structure_template_node,
#                [('out_structure_list', 'in_structure')]),
#            # link objects to structures (so that each set
#            # {[object],struct} can be looped over)
#            (
#                brain_nii_node, link_objects_structures_node,
#                [('out_file', 'in_brain_list')]),
#            (
#                self.input_node, link_objects_structures_node,
#                [('input_segmentations', 'in_segmentation_list')]),
#            (
#                get_affine_node, link_objects_structures_node,
#                [('out_matrix_file', 'in_affine_list')]),
#            (
#                process_structure_node, link_objects_structures_node,
#                [('out_structure_list', 'in_structure_list')]),
#            # get structure IDs from structures
#            (
#                link_objects_structures_node, get_structureids_node,
#                [('out_structure_list', 'in_structure')]),
#            # mesh tesselation
#            (
#                link_objects_structures_node, tesselate_node,
#                [('out_segmentation_list', 'in_file')]),
#            (
#                get_structureids_node, tesselate_node,
#                [('out_structure_id', 'label_value')]),
#            # mesh smoothing
#            (
#                tesselate_node, smooth_node,
#                [('surface', 'in_file')]),
#            # mesh conversion to .vtk
#            (
#                smooth_node, get_vtk_surface_node,
#                [('surface', 'in_file')]),
#            # re-orienting meshes with affine transformation
#            (
#                get_vtk_surface_node, apply_affine_node,
#                [('converted', 'in_segmentation')]),
#            (
#                link_objects_structures_node, apply_affine_node,
#                [('out_affine_list', 'in_affine')]),
#            (
#                link_objects_structures_node, apply_affine_node,
#                [('out_brain_list', 'in_brain')]),
#            (
#                link_objects_structures_node, apply_affine_node,
#                [('out_structure_list', 'in_structure')]),
#            # get per-subject list of meshes
#            (
#                apply_affine_node, get_persubject_meshes_node,
#                [('out_affine_reoriented', 'in_mesh_list')]),
#            (
#                link_objects_structures_node, get_persubject_meshes_node,
#                [('out_map', 'in_map')]),
#            # rigid registration of structure to Colin 27 counterpart
#            (
#                get_persubject_meshes_node, rigid_reg_node,
#                [('out_persubject_mesh_list', 'in_persubject_mesh_list')]),
#            (
#                get_colin27_structure_template_node, rigid_reg_node,
#                [('out_colin27_structure_template', 'in_colin27_mesh')]),
#            (
#                process_structure_node, rigid_reg_node,
#                [('out_structure_list', 'in_structure_list')]),
#            # produce list of datasink base directories for
#            # {subject,structure} specific files
#            (
#                process_structure_node, link_basedir_structures_node,
#                [('out_structure_list', 'in_structure_list')]),
#            # Save {subject,structure} specific files to CAPS dir
#            (
#                link_basedir_structures_node, subject_structure_datasink,
#                [('out_basedir_list', 'base_directory')]),
#            (
#                get_vtk_surface_node, subject_structure_datasink,
#                [('converted', '{0}.@initial_surface'.format(data_subfolder))]),
#            (
#                apply_affine_node, subject_structure_datasink,
#                [('out_affine_reoriented', '{0}.@affine_surface'.format(data_subfolder))]),
#            # Save subject specific files to CAPS dir
#            (
#                get_affine_node, subject_datasink,
#                [('out_matrix_file', '{0}.@affine_matrix'.format(data_subfolder))]),
#            (
#                get_affine_node, subject_datasink,
#                [('out_file', '{0}.@affine_volume'.format(data_subfolder))]),
#            (
#                rigid_reg_node, subject_datasink,
#                [('out_persubject_regmesh_list', '{0}.@rigid_surface'.format(data_subfolder))]),
#            # Save structure specific files to CAPS dir
#            (
#                get_colin27_structure_template_node, structure_datasink,
#                [('out_colin27_structure_template', 'templates')]),
#            # Save group files (parameter files, list of structures
#            # file) to CAPS dir
#            (
#                process_structure_node, group_datasink,
#                [('out_structure_file', 'structures.@file')])
            ])
