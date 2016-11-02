function clinicasurfstat( inputdir, outputdir, csvfile, linearmodel, contrast, strformat, varargin)
% Saves all the output images for group analysis of T1 images smoothed data
%
% Usage: [some outputs] = clinicasurfstat( inputdir, outputdir, csvfile, linearmodel, contrast, strformat, varargin)
%
% - inputdir:  the output file from recon-all pipeline,specifically, files: ?h.thickness.fwhm**.mgh.
% - outputdir: the directory to contain the result images.
% - linearmodel: string, the linear model that fit into the GLM, for example '1+Lable'
% - contrast: string, the contrast that you want to use in the GLM, but the contrast should be inclued into the linearmodel, otherwise, you will get errors.
% - csvfile: string, the path to your csv file
% - strstrformat: string, the strstrformat that you want to use for your CSV file column variables, it depends on the CSV file.




% the following are optional, we use varargin to define the optional inputs!
% - sizeoffwhm: fwhm for the surface smoothing, default is 20, integer.
% - thresholduncorrectedpvalue: threshold to display the uncorrected Pvalue, float.
% - thresholdcorrectedpvalue: the threshold to display the corrected cluster, default is 0.05, float.
% - clusterthreshold: threshold to define a cluster in the process of cluster-wise correction, default is 0.001, float.
%
% Saves output images to outputdir as .jpg format.
%
% Note: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%  (c) Alexandre ROUTIER, Junhao WEN 
%  Written 16/06/2016

%% define the default value for inputs
if nargin < 6
    error('the function at least has 5 required inputs!');
end
if nargin == 7
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
addpath(strcat(fileparts(which(mfilename())), '/matlab-freesurfer-5.3'));
addpath(strcat(fileparts(which(mfilename())), '/tools'));
surfstathome = fileparts(which(mfilename()));
%% Load the data
fid = fopen(csvfile, 'r');
linevar = fgetl(fid);
firstline = regexp(linevar,'\s+','split'); % in matlab2016, regexp is deprecated, we use strsplit;

lencolumn = length(firstline);
if lencolumn <2
    error('requires at least 2 inputs')
end
if  strcmp(firstline{1},'ID') ~= 1
    error('the first colomn of CSV file should always be named by ID')
end

csvdata = textscan( fid, strformat, 'HeaderLines', 0);
fclose(fid);

%% read the thickness for all the subjects!
nrsubject = length(csvdata{1});
% nrfactor1 = 0; nrfactor2 = 0;
csvheader = firstline;
for indexsvheader = 2 : lencolumn
    csvheader{indexsvheader} = cell(1,3);
end

for indexsubject = 1 : nrsubject
    subjectname = csvdata{1}{indexsubject};
    interpath = strcat(inputdir, '/', subjectname, '/*/*/*/*/surf' )
    [surfsubdir, xuuu] = glob(interpath)
    surfsubdir = char(surfsubdir)
    Y = SurfStatReadData( { ...
        [ surfsubdir '/lh.thickness.fwhm' num2str(sizeoffwhm) '.fsaverage.mgh' ],...
        [ surfsubdir '/rh.thickness.fwhm' num2str(sizeoffwhm) '.fsaverage.mgh' ]} );
    %Y = SurfStatReadData( { ...
     %   [ inputdir '/' subjectname '/surf/lh.thickness.fwhm' num2str(sizeoffwhm) '.fsaverage.mgh' ],...
      %  [ inputdir '/' subjectname '/surf/rh.thickness.fwhm' num2str(sizeoffwhm) '.fsaverage.mgh' ]} );
    if size(Y, 1) ~= 1
        error('Unexpected dimension of Y in SurfStatReadData')
    end
    disp(['The subject ID is: ', subjectname] )
    if indexsubject == 1
        thicksubject = zeros(nrsubject, size(Y,2));
        indexunique = strfind(firstline, contrast);
        indexunique = find(not(cellfun('isempty', indexunique)));
        if iscell(csvdata{indexunique})
            uniquelabels = unique(csvdata{indexunique}); factor1 = char(uniquelabels(1)); factor2 = char(uniquelabels(2));
            if length(uniquelabels) ~= 2
                error('there should be just 2 different groups!')
            end
        end
    end
    thicksubject(indexsubject, :) = Y;
         
    % create the csvLabel to display the Population info, but if it is
    % necessary to do it here?
%     if  strcmp(csvdata{indexunique}{indexsubject},factor1)
%         for indexlabel = 2 : (lencolumn)
%             csvheader{indexlabel}{1}(end+1)  = csvdata{indexlabel}(indexsubject);
%         end
%         nrfactor1 = nrfactor1 + 1;
%     elseif strcmp(csvdata{indexunique}{indexsubject}, factor2)
%         for indexlabel = 2 : (lencolumn)
%             csvheader{indexlabel}{1}(end+1)  = csvdata{indexlabel}(indexsubject);
%         end
%         nrfactor2   = nrfactor2 + 1;
%     else
%         error(['uniquelabels should be ' factor1 ' or ' factor2 '.'])
%     end
end

% Sort our label and all the other columns
% csvsorted = csvdata;
% j = 1; k = 1;
% for indexsorted= 1 : nrsubject
%     if strcmp(csvdata{indexunique}{indexsorted},factor1)
%         for index = 1 : lencolumn
%             csvsorted{index}(j) = csvdata{index}(indexsorted);
%         end
%         j = j+1;
%     elseif strcmp(csvdata{indexunique}{indexsorted},factor2)
%         for index = 1 : lencolumn
%             csvsorted{index}(k + nrfactor1) = csvdata{index}(indexsorted);
%         end
%         k = k+1;
%     else
%         error(['uniquelabels should be ' factor1 'or' factor2 ').'])
%     end
% end
% 
% %% Population info, here, we can have a if condition and create a function to cal the population info
% disp('###')
% disp(['The number of subjectes is: ' num2str(nrsubject) ])
% disp(['The number of ' factor1 ' is: ' num2str(nrfactor1) ])
% disp(['The number of ' factor2 ' is: ' num2str(nrfactor2) ])
% 
% for indexpop = 1 : lencolumn
%     if  iscell(csvsorted{indexpop}) == 0 % the 2th cell contains the continue factor popInfo
%         csvheader{indexpop}{2}(end+ 1) =  mean(csvsorted{indexpop}(1: nrfactor1));
%         disp([factor1 ' : the mean of ' num2str(indexpop) 'th facor is ' num2str(csvheader{indexpop}{2}(end))])
%         csvheader{indexpop}{2}(end+ 1) =  mean(csvsorted{indexpop}(nrfactor1 +1 : nrsubject));
%         disp([ factor2 ' : the mean of ' num2str(indexpop) 'th facor is ' num2str(csvheader{indexpop}{2}(end))])
%         csvheader{indexpop}{3}(end+ 1) =  std(csvsorted{indexpop}(1: nrfactor1));
%         disp([factor1 ' : the sd of ' num2str(indexpop) 'th facor is ' num2str(csvheader{indexpop}{3}(end))])
%         csvheader{indexpop}{3}(end+ 1) =  std(csvsorted{indexpop}(nrfactor1 +1 : nrsubject));
%         disp([ factor2 ' : the sd of ' num2str(indexpop) 'th facor is ' num2str(csvheader{indexpop}{3}(end))])
%         
%     end
% end
% disp('###')

%% Load average surface & creation of the mask :
averagesurface = SurfStatReadSurf( { strcat(surfstathome,'/fsaverage/lh.pial') , strcat(surfstathome,'/fsaverage/rh.pial') } );
thicksubject = thicksubject';
mask = thicksubject(:,1)>0;
mask = mask';

% create the Term that will be defined as contrast
for i = 1:length(csvdata)-1
    eval([firstline{i+1} '= term(csvdata{i+1});'])
end

%% Create Images folder for output
cd(outputdir)
if exist('ClinicaSurfStatOutput', 'dir') ~= 7
    mkdir ClinicaSurfStatOutput
end

%% Convert the data into SurfStat

if iscell(csvdata{indexunique})
    factorname = char(eval(contrast));
    if factor1 == factorname{1}
        contrastpos    = eval([contrast '(1)']) - eval([ contrast '(2)']); % use char(eval(contrast))
        contrasteffectgroupneg    = -contrastpos;
    else
        contrastpos    = eval([contrast '(2)']) - eval([ contrast '(1)']); % use char(eval(contrast))
        contrasteffectgroupneg    = -contrastpos;
    end

    thicksubject = thicksubject';

    slmmodel  = SurfStatLinMod(thicksubject, eval(linearmodel), averagesurface);
    disp(['The GLM linear model is: ', linearmodel])

    %% Clear the variables which will not be used later
    clearvars csvsorted csvdata thicksubject

    % Contrast Positive:
    tic;
    slm = SurfStatT( slmmodel, contrastpos );
    SurfStatView( slm.t .* mask, averagesurface, [ 'ContrastPo-value of the T-statistic for ' factor1 '-' factor2]);
    save2jpeg(strcat(outputdir,'/ClinicaSurfStatOutput/ContrastPositive-TValue.jpg'));
    disp('Contrast Positive: Tvalue'); toc;

    % Computation of the uncorrected p-value:
    tic;
    %Method 1:
    %t = tpdf(abs(slm.t), slm.df);
    %uncorrectedpValues = double(t<=0.05).*t+double(t>0.05);
    %Method 2:
    uncorrectedpValues = 2*(1-tcdf(abs(slm.t),slm.df)); % here, we consider it to be 2-tailed t-distribution
    clearvars struct; struct.P = uncorrectedpValues; struct.mask = mask; struct.thresh = thresholduncorrectedpvalue;
    SurfStatView( struct, averagesurface, [ '(ContrastPo-Uncorrected P-values)' factor1 '-' factor2 ]);
    save2jpeg(strcat(outputdir, '/ClinicaSurfStatOutput/ContrastPositive-UncorrectedPValue.jpg'));   
    disp('Contrast Positive: uncorrected Pvalue'); toc;

    % Computation of the corrected p-values: P-value threshold or statistic threshold for defining clusters, 0.001 by default
    tic;
    [ pval, ~, clus ] = SurfStatP( slm , mask, clusterthreshold);
    pval.thresh = thresholdcorrectedpvalue;
    SurfStatView( pval, averagesurface, ['(ContrastPo-Corrected P-values) ' factor1 '-' factor2 ' (clusterthreshold = ' num2str(clusterthreshold) ')']);
    save2jpeg(strcat('ClinicaSurfStatOutput/ContrastPositive-CorrectedPValue.jpg' ));
    disp('Contrast Positive: Corrected Pvalue'); toc;

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
    save2jpeg(strcat('ClinicaSurfStatOutput/ContrastPositive-FalseDiscoveryRate.jpg')) ;
    disp('Contrast Positive: FDR'); toc;
    
    %% Contrast Negative:
    % Computation of the T-statistic̣:  T statistics for a contrast in a univariate or multivariate model.
    tic;
    slm = SurfStatT( slmmodel, contrasteffectgroupneg );
    SurfStatView( slm.t .* mask, averagesurface, [ 'ContrastNe-value of the T-statistic for ' factor1 '-' factor2 ]);
    save2jpeg( strcat('ClinicaSurfStatOutput/ContrastNegative-TValue.jpg'));
    disp('Contrast Negative: Tvalue'); toc;
    
    % Computation of the uncorrected p-value:
    tic;
    uncorrectedpValues = 2*(1-tcdf(abs(slm.t),slm.df)); % here, we consider it to be 2-tailed t-distribution
    clearvars struct; struct.P = uncorrectedpValues; struct.mask = mask; struct.thresh = thresholduncorrectedpvalue;
    SurfStatView( struct, averagesurface, [ '(ContrastNe-Uncorrected P-values )' factor1 '-' factor2]);
    save2jpeg(strcat('ClinicaSurfStatOutput/ContrastNegative-UncorrectedPValue.jpg'));
    disp('Contrast Negative: Uncorrected Pvalue'); toc;
    

    % Computation of the corrected p-values: P-value threshold or statistic threshold for defining clusters, 0.001 by default
    tic;
    [ pval, ~, clus ] = SurfStatP( slm , mask, clusterthreshold);
    pval.thresh = thresholdcorrectedpvalue;
    %SurfStatView( pval, averagesurface, ['(ContrastNe-Corrected P-values )' factor1 '-' factor2 ' (clusterthreshold = ' num2str(clusterthreshold) ')']);
    % to change the background color(black), uncomment this line.
    SurfStatView( pval, averagesurface, ['(ContrastNe-Corrected P-values )' factor1 '-' factor2 ' (clusterthreshold = ' num2str(clusterthreshold) ')'], 'black');

    save2jpeg(strcat('ClinicaSurfStatOutput/ContrastNegative-CorrectedPValue.jpg'));
    disp('Contrast Negative: Corrected Pvalue'); toc;

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
    save2jpeg(strcat('ClinicaSurfStatOutput/ContrastNegative-FalseDiscoveryRate.jpg' ));
    disp('Contrast Negative: FDR'); toc;
    
else
    contrastpos    = eval(contrast);
    
    thicksubject = thicksubject';

    slm  = SurfStatLinMod(thicksubject, eval(linearmodel), averagesurface);
    disp(['The GLM linear model is: ', linearmodel])

    %% Clear the variables which will not be used later
    clearvars csvsorted csvdata thicksubject

    % Contrast Positive:
    tic;
    slm = SurfStatT( slm, contrastpos );
    SurfStatView( slm.t .* mask, averagesurface, [ 'T-statistic for ' contrast ]);
    save2jpeg(strcat(outputdir,'/ClinicaSurfStatOutput/TValue.jpg'));
    disp('Tvalue'); toc;

    % Computation of the uncorrected p-value:
    tic;
    uncorrectedpValues = 2*(1-tcdf(abs(slm.t),slm.df)); % here, we consider it to be 2-tailed t-distribution
    clearvars struct; struct.P = uncorrectedpValues; struct.mask = mask; struct.thresh = thresholduncorrectedpvalue;
    SurfStatView( struct, averagesurface, [ '(Uncorrected P-values)' contrast ]);
    save2jpeg(strcat(outputdir, '/ClinicaSurfStatOutput/UncorrectedPValue.jpg'));
    disp('Uncorrected Pvalue'); toc;

    % Computation of the corrected p-values: P-value threshold or statistic threshold for defining clusters, 0.001 by default
    tic;
    [ pval, ~, clus ] = SurfStatP( slm , mask, clusterthreshold);
    pval.thresh = thresholdcorrectedpvalue;
    SurfStatView( pval, averagesurface, ['(Corrected P-values) ' contrast ' (clusterthreshold = ' num2str(clusterthreshold) ')']);
    save2jpeg(strcat('ClinicaSurfStatOutput/CorrectedPValue.jpg' ));
    disp('Corrected Pvalue'); toc;

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
    save2jpeg(strcat('ClinicaSurfStatOutput/FalseDiscoveryRate.jpg')) ;
    disp('FDR'); toc;

end

