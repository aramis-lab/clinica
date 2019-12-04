spm('defaults','pet');
spm_jobman('initcfg');

% Setup contrast task
weights = [-1 1 [0] [0]];
matlabbatch{1}.spm.stats.con.spmmat = {@SPMMAT};
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'group1_less_than_group2';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = weights;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;

weights2 = [1 -1 [0] [0]];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'group2_less_than_group1';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = weights2;
matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

spm_jobman('run', matlabbatch)


% MAKE IT 2 files !!

% Setup resut task
matlabbatch{1}.spm.stats.results.spmmat = {@SPMMAT};
matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
matlabbatch{1}.spm.stats.results.conspec.contrasts = Inf;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{1}.spm.stats.results.conspec.extent = 0;
matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.export{1}.png = true;


% RUN BATCH JOB
spm_jobman('run', matlabbatch)
