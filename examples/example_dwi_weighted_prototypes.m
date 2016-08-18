close all
clearvars
clc

cd('/Users/pietro.gori/Softwares/clinica/examples/external-data/WeightedPrototypes')

filename_bundle='../Bundle_small.vtk';
lambda_g=3;
lambda_a=3;
lambda_b=3;
path_matlab_functions='/Users/pietro.gori/Softwares/clinica/clinica/lib/WeightedPrototypes/Matlab_Functions';
path_CPP_code='/Users/pietro.gori/Softwares/clinica/clinica/lib/WeightedPrototypes/CPP_code';
path_Community_latest='/Users/pietro.gori/Softwares/clinica/clinica/lib/WeightedPrototypes/Community_latest';
bound_limit_input=[];
degree_precision_input=[];
num_iter_modularity_input=[];
minimum_number_fibers_cluster_input=1;
minValueTau_input=[];
increase_radius_input=[];

dossier=['Example_bundle_small_lambda_g_' num2str(lambda_g) '_lambda_a_' num2str(lambda_a) '_lambda_b_' num2str(lambda_b)];

if exist(dossier,'dir')
    error('Example already computed. Rename/delete the folder if you want to run it again.')
end

mkdir(dossier)
cd(dossier)

weighted_prototypes(filename_bundle,lambda_g,lambda_a,lambda_b,path_matlab_functions,path_CPP_code,path_Community_latest,bound_limit_input,degree_precision_input,num_iter_modularity_input,minimum_number_fibers_cluster_input,minValueTau_input,increase_radius_input)

