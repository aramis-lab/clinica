"""

"""


from clinica.pipeline.t1.t1_freesurfer import datagrabber_t1_freesurfer_pipeline

reconall_wf = datagrabber_t1_freesurfer_pipeline('/Users/jeremy.guillon/Repositories/multiconproject/data/HMTC_BIDS',
                                                 '/Users/jeremy.guillon/Repositories/multiconproject/data/HMTC_CAPS',
                                                 recon_all_args='-qcache',
                                                 working_directory='/Users/jeremy.guillon/Tmp3')
reconall_wf.run()
