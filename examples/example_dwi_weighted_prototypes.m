close all
clearvars
clc

cd('/Users/pietro.gori/Softwares/clinica/examples/external-data/dwi_weighted_prototypes')

addpath('/Users/pietro.gori/Softwares/clinica/clinica/lib/weighted_prototypes_lib')

filename_bundle='../Bundle_small.vtk';
lambda_g=3;
lambda_a=3;
lambda_b=3;
path_matlab_functions='/Users/pietro.gori/Softwares/clinica/clinica/lib/weighted_prototypes_lib/matlab_functions';
path_cpp_code='/Users/pietro.gori/Softwares/clinica/clinica/lib/weighted_prototypes_lib/cpp_code';
path_community_latest='/Users/pietro.gori/Softwares/clinica/clinica/lib/weighted_prototypes_lib/community_latest';
bound_limit_input=0;
degree_precision_input=0;
num_iter_modularity_input=0;
minimum_number_fibers_cluster_input=1;
minValueTau_input=0;
increase_radius_input=0;

dossier=['Example_bundle_small_lambda_g_' num2str(lambda_g) '_lambda_a_' num2str(lambda_a) '_lambda_b_' num2str(lambda_b)];

if exist(dossier,'dir')
    error('Example already computed. Rename/delete the folder if you want to run it again.')
end

mkdir(dossier)
cd(dossier)

working_dir=pwd;
weighted_prototypes(working_dir,filename_bundle,lambda_g,lambda_a,lambda_b,path_matlab_functions,path_cpp_code,path_community_latest,bound_limit_input,degree_precision_input,num_iter_modularity_input,minimum_number_fibers_cluster_input,minValueTau_input,increase_radius_input)

