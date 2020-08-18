function clinicasurfstat( inputdir, outputdir, tsvfile, designmatrix, contrast, strformat, glmtype, grouplabel, freesurferhome, surface_file, feature_label, varargin)
% Saves all the output images for group analysis of T1 images smoothed data
%
% Usage: [some outputs] = clinicasurfstat( inputdir, outputdir, tsvfile, designmatrix, contrast, strformat, glmtype, varargin)
%
% - inputdir:  the output file from recon-all pipeline,specifically, files: ?h.thickness.fwhm**.mgh.
% - outputdir: the directory to contain the result images.
% - designmatrix: string, the linear model that fit into the GLM, for example '1+Lable'
% - contrast: string, the contrast that you want to use in the GLM, but the contrast should be inclued into the designmatrix, otherwise, you will get errors.
% - tsvfile: string, the path to your tsv file
% - strstrformat: string, the strstrformat that you want to use for your tsv file column variables, it depends on the tsv file.




% the following are optional, we use varargin to define the optional inputs!
% - sizeoffwhm: fwhm for the surface smoothing, default is 20, integer.
% - thresholduncorrectedpvalue: threshold to display the uncorrected Pvalue, float.
% - thresholdcorrectedpvalue: the threshold to display the corrected cluster, default is 0.05, float.
% - clusterthreshold: threshold to define a cluster in the process of cluster-wise correction, default is 0.001, float.
%
% Saves output images to outputdir as .jpg format.
%
% Note: the default threshold for RFT(corrected pvalue) is 0.001 and for FDR is 0.05. The RFT one is extremely strict, and Boris Bernhardt and Keith Worsley have suggested that it may be loosened to 0.01, maybe up to to 0.025.

%  (c)Junhao Wen, Alexandre ROUTIER
%  Written 16/06/2016

%% define the default value for inputs
if nargin < 9
    error('the function at least has 9 required inputs!');
end
if nargin == 10
    error('Number of inputs is wrong!')
end
sizeoffwhm = 20;
thresholduncorrectedpvalue = 0.001;
thresholdcorrectedpvalue = 0.050;
clusterthreshold = 0.001;
pp = varargin;
while length(pp) >= 2,
    prop = pp{1};
    val = pp{2};
    pp = pp(3 : end);
    if ischar(prop)
        switch prop
            case 'sizeoffwhm'
                sizeoffwhm = val;
            case 'thresholduncorrectedpvalue'
                thresholduncorrectedpvalue = val;
            case 'thresholdcorrectedpvalue'
                thresholdcorrectedpvalue = val;
            case 'clusterthreshold'
                clusterthreshold = val;
            otherwise
                error('Optional inputs are wrong, we do not have this optional input!')
        end
    else
        error('The type of default variables are not correct!')
    end
end

disp('###For OpengGl, we should choose Hardware to render the images###');
opengl info

%% define the path
addpath(fileparts(which(mfilename())));
addpath(strcat(fileparts(which(mfilename())), '/SurfStat'));
fsaveragepath = strcat(freesurferhome, '/subjects/fsaverage/surf');
addpath(strcat(freesurferhome, '/matlab'));
surfstathome = fileparts(which(mfilename()));
%% Load the data
fid = fopen(tsvfile, 'r');
linevar = fgetl(fid);
firstline = regexp(linevar,'\s+','split'); % in matlab2016, regexp is deprecated, we use strsplit;

lencolumn = length(firstline);
if lencolumn <2
    error('requires at least 2 inputs columns')
end
if  strcmp(firstline{1},'participant_id') ~= 1
    error('the first colomn of TSV file should always be named by participant_id')
end
if  strcmp(firstline{2},'session_id') ~= 1
    error('the first colomn of TSV file should always be named by session_id')
end
tsvdata = textscan( fid, strformat, 'HeaderLines', 0);
fclose(fid);

%% read the files for all the subjects!
nrsubject = length(tsvdata{1});
% nrfactor1 = 0; nrfactor2 = 0;
csvheader = firstline;
for indexsvheader = 2 : lencolumn
    csvheader{indexsvheader} = cell(1,3);
end

for indexsubject = 1 : nrsubject
    subjectname = tsvdata{1}{indexsubject};
    sessionname = tsvdata{2}{indexsubject};
    surf_file = regexprep(surface_file, '@subject', subjectname);
    surf_file = regexprep(surf_file, '@session', sessionname);
    surf_file = regexprep(surf_file, '@fwhm', num2str(sizeoffwhm));
    surf_file_lh = fullfile(inputdir, regexprep(surf_file, '@hemi', 'lh'))
    surf_file_rh = fullfile(inputdir, regexprep(surf_file, '@hemi', 'rh'))
    Y = SurfStatReadData({surf_file_lh, surf_file_rh});
    if size(Y, 1) ~= 1
        error('Unexpected dimension of Y in SurfStatReadData')
    end
    disp(['The subject ID is: ', subjectname, '_', sessionname] )
    if indexsubject == 1 % check the header format of the input tsv
        thicksubject = zeros(nrsubject, size(Y,2));
        %if startsWith(contrast, '-') this works for matlab2016b
        if strfind(contrast, '-') % check for negative correlation analysis
            abscontrast = contrast(2:end);
        else
            abscontrast = contrast;
        end
        if strfind(contrast, '*')
            with_intercation = 1;
            disp('You include interaction as covariate in you model, please carefully check the format of your tsv files')
        else  % the case not include the interaction, to check the format of the tsv
            with_intercation = 0;
            indexunique = strfind(firstline, abscontrast);
            indexunique = find(not(cellfun('isempty', indexunique)));
            if iscell(tsvdata{indexunique})
                uniquelabels = unique(tsvdata{indexunique}); % this is to find the unique levels for the group/diagnosis, should just have two level, like CN vs PT
                if length(uniquelabels) ~= 2
                    error('For group comparison, there should be just 2 different groups!')
                end
            end
        end
    end
    thicksubject(indexsubject, :) = Y;
end

%% Load average surface & creation of the mask :
averagesurface = SurfStatReadSurf( { strcat(fsaveragepath,'/lh.pial') , strcat(fsaveragepath,'/rh.pial') } );
thicksubject = thicksubject';
mask = thicksubject(:,1)>0;  % alternatively, we can use SurfStatMaskCut to extract the mask too, but this mask includes still the brain stem
mask = mask';

% create the Term that will be defined as contrast
for i = 1:length(tsvdata)-1
    eval([firstline{i+1} '= term(tsvdata{i+1});'])
end

cd(outputdir)
measurelabel = feature_label;

%% Convert the data into SurfStat

switch glmtype
    case 'group_comparison'
        if with_intercation ~= 1 % case without interaction, contrast should be negative and positive
            contrastpos    = eval([contrast '(1)']) - eval([ contrast '(2)']); % use char(eval(contrast))
            contrasteffectgroupneg    = eval([contrast '(2)']) - eval([ contrast '(1)']); % use char(eval(contrast))
            cell_factor = char(term(tsvdata{indexunique}));
            factor1 = cell_factor{1};
            factor2 = cell_factor{2};

            thicksubject = thicksubject';
            slmmodel = SurfStatLinMod(thicksubject, eval(designmatrix), averagesurface);
            disp(['The GLM linear model is: ', designmatrix])

            %% Clear the variables which will not be used later
            clearvars csvsorted tsvdata thicksubject

            % Contrast Positive:
            tic;
            slm = SurfStatT( slmmodel, contrastpos );
            SurfStatView( slm.t .* mask, averagesurface, [ 'ContrastPo-value of the T-statistic for ' factor1 '-' factor2]);
            set(gcf,'PaperPositionMode','auto');
            print(['group-' grouplabel '_' factor2 '-lt-' factor1 '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_TStatistics'] ,'-djpeg','-r0'); close
            disp('Contrast Positive: t_value'); toc;
            tvaluewithmask = slm.t .* mask;
            save(['group-' grouplabel '_' factor2 '-lt-' factor1 '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_TStatistics.mat'],'tvaluewithmask');

            % Computation of the uncorrected p-value:
            tic;
            %Method 1:
            %t = tpdf(slm.t, slm.df);
            %uncorrectedpvalues = double(t<=0.05).*t+double(t>0.05);
            %Method 2:
            uncorrectedpvalues = 1-tcdf(slm.t,slm.df);
            clearvars struct; struct.P = uncorrectedpvalues; struct.mask = mask; struct.thresh = thresholduncorrectedpvalue;
            SurfStatView( struct, averagesurface, [ 'ContrastPo-Uncorrected P-values(' num2str(thresholduncorrectedpvalue) ')' factor1 '-' factor2 ]);
            set(gcf,'PaperPositionMode','auto');
            print(['group-' grouplabel '_' factor2 '-lt-' factor1 '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_uncorrectedPValue'] ,'-djpeg','-r0'); close
            disp('Contrast Positive: uncorrected Pvalue'); toc;
            uncorrectedpvaluesstruct = struct;
            save(['group-' grouplabel '_' factor2 '-lt-' factor1 '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_uncorrectedPValue.mat'],'uncorrectedpvaluesstruct');

            % Computation of the corrected p-values: P-value threshold or statistic threshold for defining clusters, 0.001 by default
            tic;
            [ pval, ~, clus ] = SurfStatP( slm , mask, clusterthreshold);
            pval.thresh = thresholdcorrectedpvalue;
            SurfStatView( pval, averagesurface, ['(ContrastPo-Corrected P-values) ' factor1 '-' factor2 ' (clusterthreshold = ' num2str(clusterthreshold) ')']);
            set(gcf,'PaperPositionMode','auto');
            print(['group-' grouplabel '_' factor2 '-lt-' factor1 '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_correctedPValue'] ,'-djpeg','-r0'); close
            disp('Contrast Positive: Corrected Pvalue'); toc;
            correctedpvaluesstruct = pval;
            save(['group-' grouplabel '_' factor2 '-lt-' factor1 '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_correctedPValue.mat'],'correctedpvaluesstruct');

            disp('###')
            disp('After correction(Clusterwise Correction for Multiple Comparisons): ')
            if isempty(clus) ~= 1
                disp(['#Clusters found:                             ' num2str(  length(clus.P)        ) ])
                disp(['#Significative clusters (after correction) : ' num2str(  length(find(clus.P<=thresholdcorrectedpvalue))  ) ])
            else
                disp('No cluster found!');
            end
            disp('###');

            % Computation of the false discovery rate :
            tic;
            qval = SurfStatQ( slm , mask );
            SurfStatView( qval, averagesurface, ['ContrastPo-False discovery rate ' factor1 '-' factor2 ]);
            set(gcf,'PaperPositionMode','auto');
            print(['group-' grouplabel '_' factor2 '-lt-' factor1 '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_FDR'] ,'-djpeg','-r0'); close
            disp('Contrast Positive: FDR'); toc;
            qvaluesstruct = qval;
            save(['group-' grouplabel '_' factor2 '-lt-' factor1 '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_FDR.mat'],'qvaluesstruct');

            %% Contrast Negative:
            % Computation of the T-statisticÌ£:  T statistics for a contrast in a univariate or multivariate model.
            tic;
            slm = SurfStatT( slmmodel, contrasteffectgroupneg );
            SurfStatView( slm.t .* mask, averagesurface, [ 'ContrastNe-value of the T-statistic for ' factor1 '-' factor2 ]);
            set(gcf,'PaperPositionMode','auto');
            print(['group-' grouplabel '_' factor1 '-lt-' factor2 '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_TStatistics'] ,'-djpeg','-r0'); close
            disp('Contrast Negative: t_value'); toc;
            tvaluewithmask = slm.t .* mask;
            save(['group-' grouplabel '_' factor1 '-lt-' factor2 '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_TStatistics.mat'],'tvaluewithmask');

            % Computation of the uncorrected p-value:
            tic;
            %t = tpdf(abs(slm.t), slm.df);
            %uncorrectedpvalues = double(t<=0.05).*t+double(t>0.05);
            uncorrectedpvalues = 1-tcdf(slm.t,slm.df);
            clearvars struct; struct.P = uncorrectedpvalues; struct.mask = mask; struct.thresh = thresholduncorrectedpvalue;
            SurfStatView( struct, averagesurface, [ 'ContrastNe-Uncorrected P-values(' num2str(thresholduncorrectedpvalue) ')' factor1 '-' factor2]);
            set(gcf,'PaperPositionMode','auto');
            print(['group-' grouplabel '_' factor1 '-lt-' factor2 '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_uncorrectedPValue'] ,'-djpeg','-r0'); close
            disp('Contrast Negative: Uncorrected Pvalue'); toc;
            uncorrectedpvaluesstruct = struct;
            save(['group-' grouplabel '_' factor1 '-lt-' factor2 '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_uncorrectedPValue.mat'],'uncorrectedpvaluesstruct');

            % Computation of the corrected p-values: P-value threshold or statistic threshold for defining clusters, 0.001 by default
            tic;
            [ pval, ~, clus ] = SurfStatP( slm , mask, clusterthreshold);
            pval.thresh = thresholdcorrectedpvalue;
            SurfStatView( pval, averagesurface, ['(ContrastNe-Corrected P-values )' factor1 '-' factor2 ' (clusterthreshold = ' num2str(clusterthreshold) ')']);
            % to change the background color(black), uncomment this line.
            %SurfStatView( pval, averagesurface, ['(ContrastNe-Corrected P-values )' factor1 '-' factor2 ' (clusterthreshold = ' num2str(clusterthreshold) ')'], 'black');
            set(gcf,'PaperPositionMode','auto');
            print(['group-' grouplabel '_' factor1 '-lt-' factor2 '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_correctedPValue'] ,'-djpeg','-r0'); close
            disp('Contrast Negative: Corrected Pvalue'); toc;
            correctedpvaluesstruct = pval;
            save(['group-' grouplabel '_' factor1 '-lt-' factor2 '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_correctedPValue.mat'],'correctedpvaluesstruct');

            disp('###')
            disp('After correction(Clusterwise Correction for Multiple Comparisons): ')
            if isempty(clus) ~= 1
                disp(['#Clusters found:                             ' num2str(  length(clus.P)        ) ])
                disp(['#Significative clusters (after correction) : ' num2str(  length(find(clus.P<=thresholdcorrectedpvalue))  ) ])
            else
                disp('No cluster found!')
            end
            disp('###')

            % Computation of the false discovery rate :
            tic;
            qval = SurfStatQ( slm , mask );
            SurfStatView( qval, averagesurface, ['ContrastNe-False discovery rate ' factor1 '-' factor2 ]);
            set(gcf,'PaperPositionMode','auto');
            print(['group-' grouplabel '_' factor1 '-lt-' factor2 '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_FDR'] ,'-djpeg','-r0'); close
            disp('Contrast Negative: FDR'); toc;
            qvaluesstruct = qval;
            save(['group-' grouplabel '_' factor1 '-lt-' factor2 '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_FDR.mat'],'qvaluesstruct');
        else
            disp(['The GLM linear model is: ', designmatrix])
            disp(['The contrast here is the interaction between one continue variable and one categorical variable: ', contrast])

            strs_contrast = strsplit(contrast,'*');
            if size(eval(strs_contrast{1}), [, 2]) == 1 % test the size of the continue and categorical term
                continue_contrast = strs_contrast{1};
                categorical_contrast = strs_contrast{2};
            else
                continue_contrast = strs_contrast{2};
                categorical_contrast = strs_contrast{1} ;
            end
            categorical_term = eval(categorical_contrast);
            continue_term = eval(continue_contrast);
            contrast_interaction    = continue_term(1).*categorical_term(1) - continue_term(1).*categorical_term(2);
            thicksubject = thicksubject';
            slm = SurfStatLinMod(thicksubject, eval(designmatrix), averagesurface);
            slm = SurfStatT( slm, contrast_interaction);
            SurfStatView( slm.t .* mask, averagesurface, ['T-statistic for interaction' contrast]);
            set(gcf,'PaperPositionMode','auto');
            print(['interaction-' contrast '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_TStatistics'] ,'-djpeg','-r0'); close
            tvaluewithmask = slm.t .* mask;
            save(['interaction-' contrast '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_TStatistics.mat'],'tvaluewithmask');
            % find the maximum t score for a voxel and plot the regression
            % line to check the slopes for two groups
            biggestT = find(slm.t == max( slm.t ));
            Yseed = double( thicksubject(:,biggestT(1)) );
            SurfStatPlot( continue_term, Yseed, 1, categorical_term );
            xlabel(strcat(continue_contrast, ':', num2str(biggestT(1))))
            ylabel('Yseed')
            set(gcf,'PaperPositionMode','auto');
            print(['Highest T value vertex Yseed versus ' continue_contrast] ,'-djpeg','-r0'); close

            % F statistics are obtained by comparing nested models
            % test if the contrast has the same format like in designmatrix
            if strfind(designmatrix, contrast) ~= 0 % contrast is inside designmatrix with the same order
                designmatrix_reduced_model = strrep(strrep(designmatrix, ' ', ''), strcat('+',contrast), '');
            else % if the order is not the same in designmatrix
                contrast_inverse = strsplit(contrast, '*');
                contrast = char(strcat(contrast_inverse(2), '*', contrast_inverse(1)));
                designmatrix_reduced_model = strrep(strrep(designmatrix, ' ', ''), strcat('+',contrast), '');
            end
            slm_reduce = SurfStatLinMod( thicksubject, eval(char(designmatrix_reduced_model)), averagesurface );
            slm = SurfStatF( slm, slm_reduce );
            SurfStatView( slm.t.*mask, averagesurface, ['F-statistic for interaction ' contrast]);
            set(gcf,'PaperPositionMode','auto');
            print(['interaction-' contrast '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_FStatistics'] ,'-djpeg','-r0'); close
            fvaluewithmask = slm.t .* mask;
            save(['interaction-' contrast '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_FStatistics.mat'],'fvaluewithmask');

            % pvalue
            [ pval, ~, clus ] = SurfStatP( slm , mask, clusterthreshold);
            pval.thresh = thresholdcorrectedpvalue;
            SurfStatView( pval, averagesurface, ['Corrected Pvalue for interaction ', contrast, ' (clusterthreshold = ' num2str(clusterthreshold) ')']);
            % to change the background color(black), uncomment this line.
            %SurfStatView( pval, averagesurface, ['(ContrastNe-Corrected P-values )' factor1 '-' factor2 ' (clusterthreshold = ' num2str(clusterthreshold) ')'], 'black');
            set(gcf,'PaperPositionMode','auto');
            print(['interaction-' contrast '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_correctedPValue'] ,'-djpeg','-r0'); close
            correctedpvaluesstruct = pval;
            save(['interaction-' contrast '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_correctedPValue.mat'],'correctedpvaluesstruct');

            % display the details of cluster-wise correction
            disp('###')
            disp('After correction(Clusterwise Correction for Multiple Comparisons): ')
            if isempty(clus) ~= 1
                disp(['#Clusters found:                             ' num2str(  length(clus.P)        ) ])
                disp(['#Significative clusters (after correction) : ' num2str(  length(find(clus.P<=thresholdcorrectedpvalue))  ) ])
            else
                disp('No cluster found!')
            end
            disp('###')

            % FDR correction with 0.05 threshold, q value
            qval = SurfStatQ( slm , mask );
            SurfStatView( qval, averagesurface, ['False discovery rate for interaction ' contrast]);
            set(gcf,'PaperPositionMode','auto');
            print(['interaction-' contrast '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_FDR'] ,'-djpeg','-r0'); close
            disp('Interaction: FDR');
            qvaluesstruct = qval;
            save(['interaction-' contrast '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_FDR.mat'],'qvaluesstruct');
        end

    case 'correlation'
        %% if the contrast is continuous variable, dont use term, just use double to fit in the test!!!
        contrastpos    = tsvdata{indexunique};

        thicksubject = thicksubject';

        slm  = SurfStatLinMod(thicksubject, eval(designmatrix), averagesurface);
        disp(['The GLM linear model is: ', designmatrix])

        %% Clear the variables which will not be used later
        clearvars csvsorted thicksubject

        % Contrast Positive:
        tic;
        if strfind(contrast, '-')
                slm = SurfStatT( slm, -contrastpos );
                strs = strsplit(contrast,'-');
                contrast = char(strs(2));
                contrastsign = 'negative';
        else
                slm = SurfStatT( slm, contrastpos );
                contrastsign = 'positive';
        end
        SurfStatView( slm.t .* mask, averagesurface, [ 'T-statistic for ' contrast ]);
        set(gcf,'PaperPositionMode','auto');
        print(['group-' grouplabel '_correlation-' contrast '_contrast-' contrastsign '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_TStatistics'] ,'-djpeg','-r0'); close
        disp('t_value'); toc;
        % to save the T value map as a .mat file for visualize
        tvaluewithmask = slm.t .* mask;
        save(['group-' grouplabel '_correlation-' contrast '_contrast-' contrastsign '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_TStatistics.mat'],'tvaluewithmask');

        % Computation of the uncorrected p-value:
        tic;
        %t = tpdf(abs(slm.t), slm.df);
        %uncorrectedpvalues = double(t<=0.05).*t+double(t>0.05);
        uncorrectedpvalues = 1-tcdf(slm.t,slm.df);
        clearvars struct; struct.P = uncorrectedpvalues; struct.mask = mask; struct.thresh = thresholduncorrectedpvalue;
        SurfStatView( struct, averagesurface, [ 'Uncorrected P-values(' num2str(thresholduncorrectedpvalue) ')' contrast ]);
        set(gcf,'PaperPositionMode','auto');
        print(['group-' grouplabel '_correlation-' contrast '_contrast-' contrastsign '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_uncorrectedPValue'] ,'-djpeg','-r0'); close
        disp('Uncorrected Pvalue'); toc;
        uncorrectedpvaluesstruct = struct;
        save(['group-' grouplabel '_correlation-' contrast '_contrast-' contrastsign '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_uncorrectedPValue.mat'],'tvaluewithmask');

        % Computation of the corrected p-values: P-value threshold or statistic threshold for defining clusters, 0.001 by default
        tic;
        [ pval, ~, clus ] = SurfStatP( slm , mask, clusterthreshold);
        pval.thresh = thresholdcorrectedpvalue;
        SurfStatView( pval, averagesurface, ['(Corrected P-values) ' contrast ' (clusterthreshold = ' num2str(clusterthreshold) ')']);
        set(gcf,'PaperPositionMode','auto');
        print(['group-' grouplabel '_correlation-' contrast '_contrast-' contrastsign '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_correctedPValue'] ,'-djpeg','-r0'); close
        disp('Corrected Pvalue'); toc;
        correctedpvaluesstruct = pval;
        save(['group-' grouplabel '_correlation-' contrast '_contrast-' contrastsign '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_correctedPValue.mat'],'tvaluewithmask');

        disp('###')
        disp('After correction(Clusterwise Correction for Multiple Comparisons): ')
        if isempty(clus) ~= 1
            disp(['#Clusters found:                             ' num2str(  length(clus.P)        ) ])
            disp(['#Significative clusters (after correction) : ' num2str(  length(find(clus.P<=thresholdcorrectedpvalue))  ) ])
        else
            disp('No cluster found!');
        end
        disp('###');

        % Computation of the false discovery rate :
        tic;
        qval = SurfStatQ( slm , mask );
        SurfStatView( qval, averagesurface, ['False discovery rate ' contrast ]);
        set(gcf,'PaperPositionMode','auto');
        print(['group-' grouplabel '_correlation-' contrast '_contrast-' contrastsign '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_FDR'] ,'-djpeg','-r0'); close
        disp('FDR'); toc;
        qvaluesstruct = qval;
        save(['group-' grouplabel '_correlation-' contrast '_contrast-' contrastsign '_measure-' measurelabel '_fwhm-' num2str(sizeoffwhm) '_FDR.mat'],'tvaluewithmask');

    otherwise
        error('Check out if you define the glmtype flag correctly, or define your own general linear model, e,g MGLM')
end
