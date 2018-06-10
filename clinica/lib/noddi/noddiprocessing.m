function noddiprocessing(output_dir, noddi_img, brain_mask, roi_mask, bval, bvec, prefix, bStep, noddi_toolbox_dir, num_cores)
%This is a function to fit NODDI onto the multiple shells diffusion data
%based on this paper: NODDI: Practical in vivo neurite orientation
%dispersion and density imaging of the human brain, Neuroimage, 2012, Gary
%Zhang.
%output_dir is the path to store the result of fitting Noddi model.
%noddi_img is the path to the preprocessing noddi nii.
%brain_mask is the path to the mask which include just the voxels inside
%the brain.
%roi_mask is the path to the a maks to include just the voxels inside a
%ROI, normally, if you care about all the brain voxels, this should be
%equal to the brain_mask.
%bval is the path to the bval file 
%bvec is the path to the rotated bvec in the preprocessing step
%prefix is the subject_id in BIDS version
%bStep is the rounded bval that you want for each shell in your data.
%num_cores, by default is 1 if not given.

% Note: this pipeline works only for non-compressed image, which means
% nii.gz does not fit in this pipeline, please convert the compressed image
% before fitting into this pipeline.


if (nargin<11)
    num_cores = 1;
end

% make sure the num_cores are integer 
if ischar(num_cores)
    num_cores = str2double(num_cores);
end

%% add path to the python.m file
addpath(fileparts(which(mfilename())));
addpath(noddi_toolbox_dir);

%%
cd(output_dir);
%%Convert the raw DWI volume into the required format with the function CreateROI
CreateROI(noddi_img, roi_mask, strcat(prefix, '_fitted_original.mat'));

%% Round the bvals with amico fsl2scheme func
python(strcat(fileparts(which(mfilename())),'/roundbval.py'), bval, strcat(prefix,'_rounded.bval'), bStep) 

%% Convert the FSL bval/bvec files into the required format with the function FSL2Protocol
protocol = FSL2Protocol(strcat(prefix,'_rounded.bval'), bvec);
% save the protocol for plotting use
save(strcat(prefix,'_protocol.mat'), 'protocol');
%For Camino users, the Scheme file can be converted into the required format with the function SchemeToProtocol:
%protocol = SchemeToProtocol('NODDI_protocol.scheme');

%% Create the NODDI model with MakeModel function
noddi = MakeModel('WatsonSHStickTortIsoV_B0');
save(strcat(prefix,'_model.mat'), 'noddi');

%% Fitting the NODDI model
% fit the model with one core or with multiple cores, if with multiple
% cores, need to set up the parallel preference locally in Matlab

tic;
if num_cores == 1
    batch_fitting_single(strcat(prefix, '_fitted_original.mat'), protocol, noddi, strcat(prefix, '_fitted_params.mat'));
else
    batch_fitting(strcat(prefix, '_fitted_original.mat'), protocol, noddi, strcat(prefix, '_fitted_params.mat'), num_cores); % with python wrapper, it seems to be a bug run parpool
end
t = toc;
fprintf('Fitting NODDI model on this subject takes %f hours ...\n', t./3600);

%% Convert the estimated NODDI parameters into volumetric parameter maps
SaveParamsAsNIfTI(strcat(prefix, '_fitted_params.mat'), strcat(prefix, '_fitted_original.mat'), brain_mask, prefix)


