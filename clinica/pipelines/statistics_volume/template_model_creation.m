% Obtained using SPM.
% Launch SPM menu > Stats > Factorial Design Specifications > Select two sample t-test
% Then View > Show m-code

spm('defaults','pet');
spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {@OUTPUTDIR};
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = {@SCANS1}';
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = {@SCANS2}';
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
