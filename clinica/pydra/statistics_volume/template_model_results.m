% This file has been obtained with SPM
% Launch SPM menu > Stats > Results
% Then View > Show m-code

spm('defaults','pet');
spm_jobman('initcfg');

% Setup resut task
matlabbatch{1}.spm.stats.results.spmmat = {@SPMMAT};
matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
matlabbatch{1}.spm.stats.results.conspec.contrasts = Inf;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = @CORRECTIONMETHOD;
matlabbatch{1}.spm.stats.results.conspec.thresh = @THRESH;
matlabbatch{1}.spm.stats.results.conspec.extent = 0;
matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.export{1}.png = true;

% RUN BATCH JOB
spm_jobman('run', matlabbatch)
