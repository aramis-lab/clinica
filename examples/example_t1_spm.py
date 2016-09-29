#
# def launch_XNATSource_example(server, user, password, experiment_dir, output_dir):
#     """
#     Example of use of the T1 preprocessing workflow with XNATSource as input.
#
#     :param server: XNAT server address
#     :param user: XNAT username
#     :param password: XNAT password
#     :param experiment_dir: Directory to run the workflow.
#     :param output_dir: Directory to save the resulting images of segmentation and registration processes.
#     :return:
#     """
#
#     # XNATSource recommended initialization is from a configuration file
#     # stored under restricted access to keep the database access identifiers safe
#     # This is just a dummy example.
#     # For more information on XNATSource:
#     # http://nipy.org/nipype/interfaces/generated/nipype.interfaces.io.html#xnatsource
#     xnat_input = pe.Node(nio.XNATSource(server=server,
#                                         user=user,
#                                         pwd=password,
#                                         infields=['project', 'subject']),
#                          name ='xnat_input')
#
#     # A more query can be built using the infields and the query template
#     xnat_input.inputs.query_template = '/projects/%s/subjects/%s/experiments/*/scans/*/resources/*/files'
#
#     xnat_input.inputs.project = 'test1'
#     xnat_input.inputs.subject = ['005_S_0221', '007_S_0068']
#
#     # Creating T1 preprocessing workflow
#     T1_SPM_prep_wf = t1_spm_prep_pipeline(experiment_dir, output_dir)
#
#     # Creating our own workflow
#     preproc_wf = pe.Workflow(name='preproc_wf')
#     preproc_wf.base_dir = experiment_dir
#     # Connecting the XNAT source as input to the preprocessing workflow
#     preproc_wf.connect([
#         (xnat_input, T1_SPM_prep_wf, [('outfiles', 'segmentation_wf.new_segment.channel_files')])
#     ])
#
#     # Running our workflow
#     preproc_wf.run('MultiProc', plugin_args={'n_procs': 10})
#
#
# # Experiment folders
# #data_dir = '/data/FAST_DRIVE2/samper/xnat_download'
# #experiment_dir = '/data/FAST_DRIVE2/samper/clinica/run'
# #output_dir = '/data/FAST_DRIVE2/samper/clinica/run/output'
#
#
# # To test XNATSource
# # server = 'http://134.157.198.180:8080/xnat'
# # user = 'user'
# # password = 'pwd'
#
# launch_DataGrabber_example(args.data_dir, args.experiment_dir, args.output_dir)
#
# #launch_XNATSource_example(server, user, password, experiment_dir, output_dir)

clinica run t1-spm-full-prep /Users/jorge.samper/Workspace/test/full_pipeline_data /Users/jorge.samper/Workspace/test/full_pipeline_output -wd /Users/jorge.samper/Workspace/test/full_pipeline_temp -np 10 -ti 1 2 3 4 -dt 1 2 -swu -swm -wdf True True -fwhm 12 12 12 -m False -vs 3 3 3

clinica run t1-spm-segment /Users/jorge.samper/Workspace/test/segment_data /Users/jorge.samper/Workspace/test/segment_output -wd /Users/jorge.samper/Workspace/test/segment_temp -np 10 -ti 1 2 3 4 -dt 1 2 -swu -swm -wdf True True


