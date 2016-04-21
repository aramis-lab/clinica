#!/usr/bin/python

#Import necessary modules from nipype
import os
from os.path import join as opj
import nipype.pipeline.engine as pe 
import nipype.interfaces.io as nio 
from nipype.interfaces.freesurfer.preprocess import ReconAll
 
#define the variables, like the inputs and outputs files
# Subjects

# subject_list = ['sub%03d'% i for i in range(1,2)] # the input's subject_id files
experiment_dir = '/aramis/home/wen/HAO_lab/Nipype/Data_AD' # the tutorial file to contain all the files
data_dir = opj(experiment_dir, 'xnat_download')
output_dir = opj(experiment_dir, 'recon-all_all_qcache')
subject_list = []
for dirpath, dirnames, filenames in os.walk(data_dir):
    subject_list = filenames
    break
#subject_list = subject_list[:2]
#Creat the pipelines that run the recon-all command
wf = pe.Workflow(name='reconall_workflow')
# wf.base_dir = opj(experiment_dir, 'working_dir')
wf.base_dir = os.path.abspath('/aramis/home/wen/HAO_lab/Nipype/Data_AD/working_dir')

#Grab the data
datasource = pe.Node(interface = nio.DataGrabber(infields=['subject_id'], outfields=['out_files']), name="datasource")
datasource.inputs.base_directory = data_dir 
datasource.inputs.template = '%s'
datasource.inputs.subject_id = subject_list
datasource.inputs.sort_filelist = True

#make the outputs files
if not os.path.exists(output_dir):
	os.mkdir(output_dir)

#Creat the MapNode
data_dir_out = ['sub%03d'% i for i in range(1,61)]
recon_all = pe.MapNode(interface=ReconAll(),name='recon_all', iterfield=['subject_id', 'T1_files'])
recon_all.inputs.subject_id = data_dir_out
recon_all.inputs.subjects_dir = output_dir
#recon_all.inputs.T1_files = subject_list
recon_all.inputs.directive = 'all'
recon_all.inputs.args = '-qcache'

#Connect the Node
wf.connect(datasource,'out_files', recon_all,'T1_files')
#wf.connect([(datasource, recon_all, [('out_file','T1_files')])])

#Run the wf
wf.run("MultiProc", plugin_args={'n_procs':4})

#Delete the tmp files
#os.system('rm -rf %s' %wf.base_dir)##this is system dependent
#os.rmdir(wf.base_dir)
