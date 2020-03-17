% This file has been obtained with SPM
% Launch SPM menu > Stats > Contrast Manager
% Then View > Show m-code

spm('defaults','pet');
spm_jobman('initcfg');

% Setup contrast task
weights = [-1 1 @COVARNUMBER];
matlabbatch{1}.spm.stats.con.spmmat = {@SPMMAT};
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = ['Hypothesis: ', @GROUP2, ' > ', @GROUP1];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = weights;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;

weights2 = [1 -1 @COVARNUMBER];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = ['Hypothesis: ', @GROUP1, ' > ', @GROUP2];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = weights2;
matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

spm_jobman('run', matlabbatch)
