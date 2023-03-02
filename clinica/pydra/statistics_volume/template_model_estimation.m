% Obtained using SPM.
% Launch SPM menu > Stats > Model Estimation
% Then View > Show m-code

spm('defaults','pet');
spm_jobman('initcfg');

matlabbatch{1}.spm.stats.fmri_est.spmmat = {@SPMMAT};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run', matlabbatch)