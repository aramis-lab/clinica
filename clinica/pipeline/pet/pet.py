import nipype.pipeline.engine as pe
import nipype.interfaces.utility as niu
import pandas as pd


def bids_caps_pet_pipeline(bids_directory,
                           caps_directory,
                           subjects_visits_tsv,
                           group_id,
                           analysis_series_id='default',
                           pet_type='FDG',
                           working_directory=None,
                           pvc=False,
                           fwhm_x=6,
                           fwhm_y=6,
                           fwhm_z=6):

    subjects_visits = pd.io.parsers.read_csv(subjects_visits_tsv, sep='\t')
    if list(subjects_visits.columns.values) != ['participant_id', 'session_id']:
        raise Exception('Subjects and visits file is not in the correct format.')
    subjects = list(subjects_visits.participant_id)
    sessions = list(subjects_visits.session_id)

    def create_pet_wf(bids_directory,
                      caps_directory,
                      group_id,
                      analysis_series_id,
                      pet_type,
                      working_directory,
                      pvc,
                      fwhm_x,
                      fwhm_y,
                      fwhm_z,
                      index,
                      subject_session):

        from clinica.pipeline.pet.pet_utils import create_datagrabber
        import os
        import os.path as op
        import nipype.pipeline.engine as pe
        import nipype.interfaces.io as nio
        from clinica.pipeline.pet.pet_workflows import pet_pipeline

        subject = subject_session[index]['subject']
        session = subject_session[index]['session']

        # DataGrabbers
        pet_template = '%s/%s/pet/%s_%s_task-rest_acq-' + pet_type + '_pet.nii*'
        pet_files = create_datagrabber('pet_files', [subject], [session], bids_directory, pet_template)

        t1_template = '%s/%s/anat/%s_%s_T1w.nii*'
        t1_files = create_datagrabber('t1_files', [subject], [session], bids_directory, t1_template)

        ff_base_dir = op.join(caps_directory, 'analysis-series-' + analysis_series_id, 'subjects')
        ff_template = '%s/%s/t1/spm/dartel/group-' + group_id + '/flow_fields/u_rc1%s_%s_T1w_Template.nii*'
        flow_fields = create_datagrabber('flow_fields', [subject], [session], ff_base_dir, ff_template)

        dartel_template = pe.Node(nio.DataGrabber(outfields=['out_files']), name='dartel_template')
        dartel_template.inputs.base_directory = caps_directory
        dartel_template.inputs.template = op.join('analysis-series-' + analysis_series_id, 'group-' + group_id,
                                                  't1/spm/final_template/Template_6.nii')
        dartel_template.inputs.sort_filelist = False

        # TODO put reference masks somewhere in clinica
        reference_mask = pe.Node(nio.DataGrabber(outfields=['out_files']), name='reference_mask')
        clinica_home = os.getenv("CLINICA_HOME")
        base_directory = op.join(clinica_home, 'clinica', 'resources')
        reference_mask.inputs.base_directory = base_directory
        reference_mask.inputs.sort_filelist = False
        # TODO DIFFERENT PET TYPES TO PROCESS
        if pet_type == 'FDG':
            reference_mask.inputs.template = 'mask_pons_eroded_6mm.nii'
        elif pet_type == 'AV45':
            reference_mask.inputs.template = 'mask_cerebellum+pons_eroded_6mm.nii'
        else:
            raise NotImplementedError

        suvr_datasink = pe.Node(nio.DataSink(), name='suvr_datasink')
        suvr_datasink.inputs.parameterization = False
        suvr_datasink.inputs.base_directory = op.join(caps_directory, 'analysis-series-' + analysis_series_id + '/subjects/' + subject + '/' + session + '/pet/preprocessing')

        pet_wf = pet_pipeline(working_directory=op.join(working_directory, 'individual_pipelines', subject + '-' + session), pvc=pvc)
        inputnode = pet_wf.get_node('inputnode')
        outputnode = pet_wf.get_node('outputnode')

        pet_wf.connect([(pet_files, inputnode, [('out_files', 'pet_image')]),
                        (t1_files, inputnode, [('out_files', 't1_image_native')]),
                        (flow_fields, inputnode, [('out_files', 'flow_fields')]),
                        (dartel_template, inputnode, [('out_files', 'dartel_template')]),
                        (reference_mask, inputnode, [('out_files', 'reference_mask')]),
                        (outputnode, suvr_datasink, [('suvr_pet_path', 'suvr')])]) # TODO Add a rename function

        if pvc:
            c_base_dir = op.join(caps_directory, 'analysis-series-' + analysis_series_id, 'subjects')
            c1_template = '%s/%s/t1/spm/segmentation/native_space/c1%s_%s_T1w_Template.nii*'
            c1_files = create_datagrabber('c1_files', [subject], [session], c_base_dir, c1_template)

            c2_template = '%s/%s/t1/spm/segmentation/native_space/c2%s_%s_T1w_Template.nii*'
            c2_files = create_datagrabber('c2_files', [subject], [session], c_base_dir, c2_template)

            c3_template = '%s/%s/t1/spm/segmentation/native_space/c3%s_%s_T1w_Template.nii*'
            c3_files = create_datagrabber('c3_files', [subject], [session], c_base_dir, c3_template)

            pet_wf.connect([(c1_files, inputnode, [('out_files', 'gm_image')]),
                            (c2_files, inputnode, [('out_files', 'wm_image')]),
                            (c3_files, inputnode, [('out_files', 'csf_image')])])

            inputnode.inputs.fwhm_x = fwhm_x
            inputnode.inputs.fwhm_y = fwhm_y
            inputnode.inputs.fwhm_z = fwhm_z

        pet_wf.run()
        return outputnode.get_output('suvr_pet_path')

    infosource = pe.Node(interface=niu.IdentityInterface(fields=['bids_directory',
                                                                 'caps_directory',
                                                                 'group_id',
                                                                 'analysis_series_id',
                                                                 'pet_type',
                                                                 'working_directory',
                                                                 'pvc',
                                                                 'fwhm_x',
                                                                 'fwhm_y',
                                                                 'fwhm_z',
                                                                 'index',
                                                                 'subject_session']),
                         name="infosource")

    infosource.inputs.bids_directory = bids_directory
    infosource.inputs.caps_directory = caps_directory
    infosource.inputs.group_id = group_id
    infosource.inputs.analysis_series_id = analysis_series_id
    infosource.inputs.pet_type = pet_type
    infosource.inputs.working_directory = working_directory
    infosource.inputs.pvc = pvc
    infosource.inputs.fwhm_x = fwhm_x
    infosource.inputs.fwhm_y = fwhm_y
    infosource.inputs.fwhm_z = fwhm_z

    subject_session_pairs = [{'subject': subjects[i], 'session': sessions[i]} for i in range(len(subjects))]
    infosource.inputs.subject_session = subject_session_pairs
    infosource.iterables = [('index', range(len(subject_session_pairs)))]

    create_wf = pe.Node(niu.Function(input_names=['bids_directory',
                                                  'caps_directory',
                                                  'group_id',
                                                  'analysis_series_id',
                                                  'pet_type',
                                                  'working_directory',
                                                  'pvc',
                                                  'fwhm_x',
                                                  'fwhm_y',
                                                  'fwhm_z',
                                                  'index',
                                                  'subject_session'],
                                     output_names=['output'],
                                     function=create_pet_wf),
                        name='create_wf')

    iter_pet_wf = pe.Workflow(name='iter_pet_wf')
    if working_directory is not None:
        iter_pet_wf.base_dir = working_directory

    iter_pet_wf.connect([(infosource, create_wf, [('bids_directory', 'bids_directory'),
                                                  ('caps_directory', 'caps_directory'),
                                                  ('group_id', 'group_id'),
                                                  ('analysis_series_id', 'analysis_series_id'),
                                                  ('pet_type', 'pet_type'),
                                                  ('working_directory', 'working_directory'),
                                                  ('pvc', 'pvc'),
                                                  ('fwhm_x', 'fwhm_x'),
                                                  ('fwhm_y', 'fwhm_y'),
                                                  ('fwhm_z', 'fwhm_z'),
                                                  ('index', 'index'),
                                                  ('subject_session', 'subject_session')])
                         ])
    return iter_pet_wf


