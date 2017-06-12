
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as niu



def statistic_T1(directory_input_image,participants_sessions_tsv, group_ID):
    import nibabel as nib
    import numpy as np
    import pandas
    import os.path
    import sys
    from atlas_average_utils import reading_image, read_csv_label, num_columns_csv, statistics

    directory_atlas=os.path.join(os.path.expandvars('$CLINICA_HOME'), 'clinica', 'resources', 'atlases_SPM')
    atlas_list=os.path.join(directory_atlas,'atlas_list.tsv')
    subjects_visits = pandas.io.parsers.read_csv((participants_sessions_tsv), sep='\t')
    atlas_name= pandas.io.parsers.read_csv((atlas_list), sep='\t')
    atlas_name=list(atlas_name.Name)
    subjects = list(subjects_visits.participant_id)
    sessions = list(subjects_visits.session_id)

    for i in xrange(len(subjects)):
        participant_id = subjects[i]
        session_id = sessions[i]
        image = reading_image(os.path.join(directory_input_image, 'subjects',
                                          participant_id, session_id, 't1', 'spm', 'dartel', 'group-'+group_ID, participant_id + '_'+session_id+'_T1w'+
                                          'segm-graymatter_space-..._modulated-on_probability.nii.gz'))
        if not os.path.isfile(image):
            image=reading_image(os.path.join(directory_input_image, 'analysis-series-default', 'subjects',
                                              participant_id, 'ses-bl', 't1', 'spm', 'dartel', 'group-'+group_ID, participant_id + '_'+session_id+'_T1w'+
                                              'segm-graymatter_space-..._modulated-on_fwhm'+fwhm_value+'_probability.nii.gz')) #we need to take the fwhm


        for jj in xrange(len(atlas_name)):
            atlas_id=atlas_name[jj]

            data=statistics(atlas_id,directory_atlas,image)

            output_directory=os.path.join(os.path.join(directory_input_image,'analysis-series-default','subjects',
                                            participant_id,session_id,'t1','spm','atlas_statistics'))
            os.system('mkdir -p '+output_directory)
            #label_list=read_csv_label(directory_atlas, atlas_id)

            data.to_csv(os.path.join(output_directory, participant_id + '_' + session_id + '_space-'+atlas_id+'_map-graymatter_' + 'statistics.tsv'), sep='\t', index=False)


def statistic_PET(directory_input_image,participants_sessions_tsv, group_ID):
    import nibabel as nib
    import numpy as np
    import pandas
    import os.path
    import sys
    from atlas_average_utils import reading_image, read_csv_label, num_columns_csv, statistics

    directory_atlas=os.path.join(os.path.expandvars('$CLINICA_HOME'), 'clinica', 'resources', 'atlases_SPM')
    atlas_list=os.path.join(directory_atlas,'atlas_list.tsv')
    subjects_visits = pandas.io.parsers.read_csv((participants_sessions_tsv), sep='\t')
    atlas_name= pandas.io.parsers.read_csv((atlas_list), sep='\t')
    atlas_name=list(atlas_name.Name)
    subjects = list(subjects_visits.participant_id)
    sessions = list(subjects_visits.session_id)

    for i in xrange(len(subjects)):
        participant_id = subjects[i]
        session_id = sessions[i]
        image = reading_image(os.path.join(directory_input_image, 'analysis-series-default', 'subjects',
                                            participant_id, 'ses-bl', 'pet', 'preprocessing', 'suvr_pet', 'suvr_wr'+participant_id + '_'+session_id+'_task-rest_acq-FDG_pet.nii.gz')) #to be modified
        for jj in xrange(len(atlas_name)):
            atlas_id=atlas_name[jj]
            data=statistics(atlas_id,directory_atlas,image)

            output_directory=os.path.join(os.path.join(directory_input_image, 'analysis-series-default','subjects',
                                        participant_id,session_id,'pet','atlas_statistics'))
            os.system('mkdir -p '+output_directory)
            #label_list=read_csv_label(directory_atlas, atlas_id)
            data.to_csv(os.path.join(output_directory, participant_id + '_' + session_id + '_space-'+atlas_id+'_map-fdg' + 'statistics2.tsv'), sep='\t', index=False)


def atlas_average_PET(working_directory = None,
                                      name='atlas_average_PET_wf'):

    atlas_average_PET = pe.Node(niu.Function(input_names=['directory_input_image', 'participants_sessions_tsv', 'group_ID'],
                                output_names=['output'],
                                function=statistic_PET),
                   name='atlas_average_PET')

    wf = pe.Workflow(name=name)

    if working_directory is not None:
        wf.base_dir = working_directory

    inputnode = pe.Node(niu.IdentityInterface(fields=['directory_input_image', 'directory_atlas', 'participants_sessions_tsv', 'group_ID','atlas_list']), name='inputnode', mandatory_inputs=True)
    outputnode = pe.Node(niu.IdentityInterface(fields=['myoutfile']), name='outputnode', mandatory_inputs=True)

    wf.connect([(inputnode, atlas_average_PET, [('directory_input_image', 'directory_input_image'),
                            ('participants_sessions_tsv', 'participants_sessions_tsv'), ('group_ID', 'group_ID')]),
                (atlas_average_PET, outputnode, [('output', 'myoutfile')])])

    return wf



def atlas_average_T1(working_directory = None,
                                      name='atlas_average_T1_wf'):

    atlas_average_T1 = pe.Node(niu.Function(input_names=['directory_input_image', 'participants_sessions_tsv', 'group_ID'],
                                output_names=['output'],
                                function=statistic_T1),
                   name='atlas_average_T1')

    wf = pe.Workflow(name=name)

    if working_directory is not None:
        wf.base_dir = working_directory

    inputnode = pe.Node(niu.IdentityInterface(fields=['directory_input_image', 'participants_sessions_tsv', 'group_ID']), name='inputnode', mandatory_inputs=True)
    outputnode = pe.Node(niu.IdentityInterface(fields=['myoutfile']), name='outputnode', mandatory_inputs=True)

    wf.connect([(inputnode, atlas_average_T1, [('directory_input_image', 'directory_input_image'),
     ('participants_sessions_tsv', 'participants_sessions_tsv'), ('group_ID', 'group_ID')]),
                (atlas_average_T1, outputnode, [('output', 'myoutfile')])])

    return wf

