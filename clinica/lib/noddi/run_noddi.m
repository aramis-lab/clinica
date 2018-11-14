% % PT
% output_dir='/aramis/home/wen/test/NODDI_fit_notbash/PT/matlab';
% noddi_img='/aramis/home/wen/test/noddi_preprocessing_test/subjects/sub-PREVDEMALS0010001BE/ses-M0/noddi/preprocessing/eddy-current-corretion/sub-PREVDEMALS0010001BE_ses-M0_eddy.nii';
% brain_mask='/aramis/home/wen/test/noddi_preprocessing_test/subjects/sub-PREVDEMALS0010001BE/ses-M0/noddi/preprocessing/sub-PREVDEMALS0010001BE_ses-M0_unwarp_B0_topup_mask.nii';
% bval='/aramis/home/wen/test/noddi_preprocessing_test/subjects/sub-PREVDEMALS0010001BE/ses-M0/noddi/preprocessing/original/sub-PREVDEMALS0010001BE_ses-M0.bval';
% bvec='/aramis/home/wen/test/noddi_preprocessing_test/subjects/sub-PREVDEMALS0010001BE/ses-M0/noddi/preprocessing/original/sub-PREVDEMALS0010001BE_ses-M0.bvec';
% prefix='matlab';
% noddi_model='WatsonSHStickTortIsoV_B0';
% num_cores=4;
% bStep ='0,300,700,2200';
% noddiprocessing(output_dir,noddi_img, brain_mask, bval, bvec, prefix, noddi_model, num_cores, bStep)

%% CN test with just one slice
output_dir='/Users/junhao.wen/test/NODDI';
noddi_img='/Users/junhao.wen/Hao/Dataset/NODDI_example_dataset/PREVDEMALS_03CN_preprocessed_data/sub-PREVDEMALS0010003PB_ses-M0_eddy.nii';
brain_mask='/Users/junhao.wen/Hao/Dataset/NODDI_example_dataset/PREVDEMALS_03CN_preprocessed_data/sub-PREVDEMALS0010003PB_ses-M0_unwarp_B0_topup_mask.nii';
roi_mask='/Users/junhao.wen/Hao/Dataset/NODDI_example_dataset/PREVDEMALS_03CN_preprocessed_data/sub-PREVDEMALS0010003PB_test_roi_mask.nii';
bval='/Users/junhao.wen/Hao/Dataset/NODDI_example_dataset/PREVDEMALS_03CN_preprocessed_data/sub-PREVDEMALS0010003PB_ses-M0.bval';
bvec='/Users/junhao.wen/Hao/Dataset/NODDI_example_dataset/PREVDEMALS_03CN_preprocessed_data/sub-PREVDEMALS0010003PB_ses-M0.bvec';
prefix='sub-PREVDEMALS0010003PB_ses-M0';
num_cores=4;
bStep ='0,300,700,2200';
noddiprocessing(output_dir,noddi_img, brain_mask, roi_mask, bval, bvec, prefix, bStep, num_cores)

%%
% %% PS
% output_dir='/aramis/home/wen/test/NODDI_fit_notbash/PS/matlab';
% noddi_img='/aramis/home/wen/test/noddi_preprocessing_test/subjects/sub-PREVDEMALS0010007RC/ses-M0/noddi/preprocessing/original/sub-PREVDEMALS0010007RC_ses-M0.nii';
% brain_mask='/aramis/home/wen/test/noddi_preprocessing_test/subjects/sub-PREVDEMALS0010007RC/ses-M0/noddi/preprocessing/sub-PREVDEMALS0010007RC_ses-M0_unwarp_B0_topup_mask.nii';
% bval='/aramis/home/wen/test/noddi_preprocessing_test/subjects/sub-PREVDEMALS0010007RC/ses-M0/noddi/preprocessing/original/sub-PREVDEMALS0010007RC_ses-M0.bval';
% bvec='/aramis/home/wen/test/noddi_preprocessing_test/subjects/sub-PREVDEMALS0010007RC/ses-M0/noddi/preprocessing/original/sub-PREVDEMALS0010007RC_ses-M0.bvec';
% prefix='matlab';
% noddi_model='WatsonSHStickTortIsoV_B0';
% num_cores=4;
% bStep ='0,300,700,2200';
% noddiprocessing(output_dir,noddi_img, brain_mask, bval, bvec, prefix, noddi_model, num_cores, bStep)
% %cd(output_dir);
% %SaveParamsAsNIfTI('FittedParams.mat', 'NODDI_roi.mat', brain_mask, prefix)
% %a=load('FittedParams.mat')
% %cd('/aramis/dataARAMIS/users/junhao.wen/PhD/PREVDEMALS/Noddi/data/NODDI_matlabtoolbox_data');
% %b=load('FittedParams.mat')
% %fprint('shit')
