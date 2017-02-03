# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 15:56:18 2016

@author: jacquemont
"""
from __future__ import absolute_import
from clinica.pipeline.statistics.connection_wise_analysis_pipeline import create_connection_wise_analysis_pipeline
from clinica.pipeline.statistics.NBS_pipeline import create_network_based_statistic_pipeline
import os
from os.path import realpath,split,join
import tempfile

data_path = join(split(realpath(__file__))[0], 'external-data/DWI_launch_preproc')

list_of_connectome_1 = 'add your path'			# List of paths to the connectome of your first group of subject
list_of_connectome_2 = 'add your path'			# List of paths to the connectome of your second group of subject
test = 't test' 					# Perform a Student T test on the connectome
FDR_correction = True 					# Correct for Family Wise Error Rate
tail = 1 						# Testing the hypothesis: Connectome group 1 > Connectome group 2. If tail set to 2 the hypothesis tested is Connectome group 1 != Connectome group 2
nb_permutation = 1000 					# Number of permutation to perform for the network-based statistics
size = 'extent' 					# Ponderate module by its extent
output_prefix = 'network-based_statistics_test' 	# Prefix used to name all network-based statistics outputs
threshold = 0.01 					# Set primary threshold to 0.01
significancy = 0.05 					# Set significancy to 0.05
save_all = True 					# Save all network-based results

working_directory = tempfile.mkdtemp()
datasink_directory = tempfile.mkdtemp()

print("Working Directory -> %s" % working_directory)
print("Datasink Directory -> %s" % datasink_directory)

if tail==1:
	print("Testing the hypothesis: Connectome group 1 > Connectome group 2")
elif tail==2:
	print("Testing the hypothesis: Connectome group 1 > Connectome group 2")
else:
	raise IOError('Set tail parametter to 1 or 2.')
	

print("Running Connection Wise analysis")
connection_wise_analysis = create_connection_wise_analysis_pipeline(list_of_connectome_1, list_of_connectome_2, test, working_directory, datasink_directory, FDR_correction=FDR_correction, tail=tail)
connection_wise_analysis.run()
print("Connection Wise analysis done")

print("Running network-based statistics")
network_based_statistic_pipeline = create_network_based_statistic_pipeline(list_of_connectome_1, list_of_connectome_2, test, nb_permutation, size, output_prefix, working_directory, datasink_directory, threshold=threshold, significancy=significancy, tail=tail, save_all=save_all)
network_based_statistic_pipeline.run()
print("Network-based statistics done")

print("Working Directory -> %s" % working_directory)
print("Datasink Directory -> %s" % datasink_directory)

