function NipypeSurfStat( ContrastLinearModel, Format, CSVFilename,PATH_TO_RECON_ALL_OUTPUTS,a_required_path, result_dir, varargin)
% Saves all the output images for group analysis of T1 images smoothed data
%
% Usage: [some outputs] = groupAnalysisReconall( ContrastLinearModel,Format, CSVFilename,PATH_TO_RECON_ALL_OUTPUTS,a_required_path,varargin);
%
% - ContrastLinearModel: string, the linear model that fit into the GLM, for example '1+Lable'
% - Format: string, the format that you want to use for your CSV file column variables, it depends on the CSV file. 
% - CSVFilename: string, the path to your csv file
% - PATH_TO_RECON_ALL_OUTPUTS:  the output file from recon-all pipeline,specifically, files: ?h.thickness.fwhm**.mgh.
% - a_required_path:  this is the path to find the this matlab function.
% - result_dir: the directory to contain the result images.


% the following are optional, we use varargin to define the optional inputs!
% - fasaverage_size_fwhm: fwhm for the surface smoothing, default is 20, integer
% - thresholdUncorrectedPValue: threshold to display the uncorrected Pvalue, float
% - thresholdCorrectedPValue: the threshold to display the corrected cluster, default is 0.05, float.
% - clusterThreshold: threshold to define a cluster in the process of cluster-wise correction, default is 0.001, float
%   
% Saves output figures to .jpg format. 
% 
% Note: here, we should define the sequences of the columns of the csv
% file, by default, it should be like this: ID, Gender, Label, age. Later,
% I will revise the function to let the user define the number and type of
% the variables that they want to add into the linear model!

% Note: when you put the path for the CSV files and the reconall output,
% be careful with the space!!!!

%  (c) Junhao WEN, anbai106@hotmail.com
%  Written 16/06/2016
%  Revised 

    %% define the default value for inputs    
    if nargin < 6
        error('the function at least has 5 required inputs!'); 
    end
    if nargin == 7
        error('Number of inputs is wrong!')
    end
    fasaverage_size_fwhm = 20;
    thresholdUncorrectedPValue = 0.001;
    thresholdCorrectedPValue = 0.05;
    clusterThreshold = 0.001;  
    pp = varargin;
    while length(pp) >= 2,
        prop = pp{1};
        val = pp{2};
        pp = pp(3 : end);
        switch prop
            case 'fasaverage_size_fwhm'
                fasaverage_size_fwhm = val;
            case 'thresholdUncorrectedPValue'
                thresholdUncorrectedPValue = val;
            case 'thresholdCorrectedPValue'
                thresholdCorrectedPValue = val;
            case 'clusterThreshold'
                clusterThreshold = val;
            otherwise
                error('Optional inputs are wrong, we do not have this optional input!')
        end
    end

    
    
     
    %% define the path
    PATH_TO_SURFSTAT          = strcat(a_required_path, '/SurfStat');
    PATH_TO_FREESURFER_MATLAB = strcat(a_required_path, '/matlab-freesurfer-5.3');   
    PATH_TO_TOOLS             = strcat(a_required_path, '/tools');
    PATH_TO_ALL_USEFUL_PACKAGE= a_required_path;
   
    %% Add the paths :
    addpath(PATH_TO_SURFSTAT)
    addpath(PATH_TO_FREESURFER_MATLAB)
    addpath(PATH_TO_TOOLS)

    %% Load the data
    fid = fopen(CSVFilename, 'r');
    lineVar = fgetl(fid); 
%     firstLine = strsplit(lineVar);% this is a new version function, If
%     your matlab version is matlab2011, you will get the error
    firstLine = regexp(lineVar,'\s+','split');
    lenthColomn = length(firstLine);
    if lenthColomn <2
        error('requires at least 2 inputs') 
    end
    if  strcmp(firstLine{1},'ID') ~= 1
        error('the first colomn of CSV file should always be named by ID')
    end 
    if  strcmp(firstLine{2},'Label') ~= 1
        error('the second colomn of CSV file should always be named by Label')
    end 
    CSVdata = textscan( fid, Format, 'HeaderLines', 0);
    fclose(fid);


    %% read the thickness for all the subjects!
    ThicknessSubjects = [];
    NumberOfSubjects = length(CSVdata{1}); 
    NumberOfLable1 = 0; NumberOfLable2 = 0;
    csvLabels = firstLine;
    for indexcsvLabels = 2 : lenthColomn
         csvLabels{indexcsvLabels} = cell(1,3);
    end
    for indexSubject = 1 : NumberOfSubjects
        nameFolder = CSVdata{1}{indexSubject};
        Y = SurfStatReadData( { ...
            [ PATH_TO_RECON_ALL_OUTPUTS '/' nameFolder '/surf/lh.thickness.fwhm' num2str(fasaverage_size_fwhm) '.fsaverage.mgh' ],...
            [ PATH_TO_RECON_ALL_OUTPUTS '/' nameFolder '/surf/rh.thickness.fwhm' num2str(fasaverage_size_fwhm) '.fsaverage.mgh' ]} );
        if size(Y, 1) ~= 1
            error('Unexpected dimension of Y in SurfStatReadData')
        end        
        ThicknessSubjects(end+1, :) = Y;
        Labels = unique(CSVdata{2}); Label1 = char(Labels(1)); Label2 = char(Labels(2));
        if length(Labels) ~= 2
            error('there should be just 2 different groups!')
        end
        disp(['#The ID: ' CSVdata{1}{indexSubject}])

        % create the csvLabel to display the Population info
        if  strcmp(CSVdata{2}{indexSubject},Label1)
            for indexData = 2 : (lenthColomn)
                csvLabels{indexData}{1}(end+1)  = CSVdata{indexData}(indexSubject);  
            end
            NumberOfLable1 = NumberOfLable1 + 1;
        elseif strcmp(CSVdata{2}{indexSubject}, Label2)
            for indexData = 2 : (lenthColomn)
                csvLabels{indexData}{1}(end+1)  = CSVdata{indexData}(indexSubject);  
            end
            NumberOfLable2   = NumberOfLable2 + 1;   
        else
            error(['Labels should be ' Label1 'or' Label2 ').'])
        end
    end

    % Sort our label and all the other columns
    CSVdataSorted = CSVdata;
    j = 1; k = 1;   
    for indexSort= 1 : NumberOfSubjects
        if strcmp(CSVdata{2}{indexSort},Label1)
            for indexC = 1 : lenthColomn
                CSVdataSorted{indexC}(j) = CSVdata{indexC}(indexSort);
            end
            j = j+1;
        elseif strcmp(CSVdata{2}{indexSort},Label2)
            for indexC = 1 : lenthColomn
                CSVdataSorted{indexC}(k + NumberOfLable1) = CSVdata{indexC}(indexSort);
            end
            k = k+1;
        else
            error(['Labels should be ' Label1 'or' Label2 ').'])
        end
    end

    %% Population info, here, we can have a if condition and create a function to cal the population info
    disp('###')
    disp(['The sorted CSV file is: ' ])
    CSVdataSorted
    disp(['The number of subjectes is: ' num2str(NumberOfSubjects) ])
    disp(['The number of ' Label1 ' is: ' num2str(NumberOfLable1) ])
    disp(['The number of ' Label2 ' is: ' num2str(NumberOfLable2) ])

    for indexPop = 1 : lenthColomn
        if  iscell(CSVdataSorted{indexPop}) == 0 % the 2th cell contains the continue factor popInfo
            csvLabels{indexPop}{2}(end + 1) =  mean(CSVdataSorted{indexPop}(1: NumberOfLable1));
            disp([ Label1 ' : the mean of ' num2str(indexPop) 'th facor is ' num2str( csvLabels{indexPop}{2}(end ))])
            csvLabels{indexPop}{2}(end + 1) =  mean(CSVdataSorted{indexPop}(NumberOfLable1 +1 : NumberOfSubjects));
            disp([ Label2 ' : the mean of ' num2str(indexPop) 'th facor is ' num2str( csvLabels{indexPop}{2}(end ))])
            csvLabels{indexPop}{3}(end + 1) =  std(CSVdataSorted{indexPop}(1: NumberOfLable1));
            disp([ Label1 ' : the sd of ' num2str(indexPop) 'th facor is ' num2str( csvLabels{indexPop}{3}(end ))])
            csvLabels{indexPop}{3}(end + 1) =  std(CSVdataSorted{indexPop}(NumberOfLable1 +1 : NumberOfSubjects));
            disp([ Label2 ' : the sd of ' num2str(indexPop) 'th facor is ' num2str( csvLabels{indexPop}{3}(end ))])
            
        end
    end
    disp('###')

    %% Load average surface & creation of the mask :
    
    cd(PATH_TO_ALL_USEFUL_PACKAGE);
    averageSurface = SurfStatReadSurf( { 'MeshesFSAVERAGE/lh.pial' , 'MeshesFSAVERAGE/rh.pial' } );
    ThicknessSubjects = ThicknessSubjects';
    mask = ThicknessSubjects(:,1)>0;
    mask = mask';
    
   %% Conversion of the data into SurfStat data
    csvTerms = cell(lenthColomn-1, 1);
    for i = 1:lenthColomn-1 
        csvTerms{i} = term(CSVdata{i + 1});
    end
    
    for i = 1:length(CSVdata)-1 
        eval([firstLine{i+1} '= term(CSVdata{i+1});'])
    end  

    contrastEffectForGroupP    = eval(['Label.' Label1]) - eval([ 'Label.' Label2 ]);
    contrastEffectForGroupN    = eval(['Label.' Label2]) - eval([ 'Label.' Label1 ]);
  
    ThicknessSubjects = ThicknessSubjects';


    slm_EffectsOnGroup  = SurfStatLinMod(ThicknessSubjects, eval(ContrastLinearModel), averageSurface);
    disp(['The GLM linear model is: ' ])
    eval(ContrastLinearModel)
    
    %% Create Images folder for output
    cd ../..
    cd examples/ClinicaSurfstat
    if exist('Figures', 'dir') ~= 7
        mkdir Figures
    end 
    cd Figures
    fileName = 'SurfStatGroupAnalysis';
    if exist(fileName,'dir') ~=7
        mkdir(fileName)
    end
    cd(PATH_TO_ALL_USEFUL_PACKAGE);

    %% Clear the variables which will not be used later
    clearvars CSVdataSorted CSVdata ThicknessSubjects
     
    % Contrast Positive: 
    slm_EffectsOnGroup = SurfStatT( slm_EffectsOnGroup, contrastEffectForGroupP );
    SurfStatView( slm_EffectsOnGroup.t .* mask, averageSurface, [ 'ContrastPo-value of the T-statistic for ' Label1 '-' Label2]);    
    cd(result_dir); 
    save2jpeg(strcat('Figures/', fileName, '/ContrastPositive-TValue.jpg')); 
    cd(PATH_TO_ALL_USEFUL_PACKAGE);

    % Computation of the uncorrected p-value:
    t = tpdf(abs(slm_EffectsOnGroup.t), slm_EffectsOnGroup.df);
    uncorrected_pValues = double(t<=0.05).*t+double(t>0.05);
    clearvars struct; struct.P = uncorrected_pValues; struct.mask = mask; struct.thresh = thresholdUncorrectedPValue;
    SurfStatView( struct, averageSurface, [ '(ContrastPo-Uncorrected P-values)' Label1 '-' Label2 ]);
    cd(result_dir); 
    save2jpeg(strcat('Figures/', fileName, '/ContrastPositive-UncorrectedPValue.jpg'));
    cd(PATH_TO_ALL_USEFUL_PACKAGE);
  
    % Computation of the corrected p-values: P-value threshold or statistic threshold for defining clusters, 0.001 by default
    [ pval, peak, clus ] = SurfStatP( slm_EffectsOnGroup , mask, clusterThreshold);
    pval.thresh = thresholdCorrectedPValue;
    SurfStatView( pval, averageSurface, ['(ContrastPo-Corrected P-values) ' Label1 '-' Label2 ' (clusterThreshold = ' num2str(clusterThreshold) ')']);
    cd(result_dir); 
    save2jpeg(strcat('Figures/', fileName, '/ContrastPositive-CorrectedPValue.jpg' )); 
    cd(PATH_TO_ALL_USEFUL_PACKAGE);
    
    disp('###')
    disp('After correction(Clusterwise Correction for Multiple Comparisons): ')
    if isempty(clus) ~= 1
        disp(['#Clusters found:                             ' num2str(  length(clus.P)        ) ])
        disp(['#Significative clusters (after correction) : ' num2str(  length(find(clus.P<=thresholdCorrectedPValue))  ) ])
    else
        disp(['No cluster found!'])
    end
    disp('###')

    % Computation of the false discovery rate :
    qval = SurfStatQ( slm_EffectsOnGroup , mask );
    SurfStatView( qval, averageSurface, ['ContrastPo-False discovery rate ' Label1 '-' Label2 ]);
    cd(result_dir); 
    save2jpeg(strcat('Figures/', fileName, '/ContrastPositive-FalseDiscoveryRate.jpg')) ;
    cd(PATH_TO_ALL_USEFUL_PACKAGE);


    %% Contrast Negative: 
    % Computation of the T-statistic̣:  T statistics for a contrast in a univariate or multivariate model.
    slm_EffectsOnGroup = SurfStatT( slm_EffectsOnGroup, contrastEffectForGroupN );    
    SurfStatView( slm_EffectsOnGroup.t .* mask, averageSurface, [ 'ContrastNe-value of the T-statistic for ' Label1 '-' Label2 ]);
    cd(result_dir); 
    save2jpeg( strcat('Figures/', fileName, '/ContrastNegative-TValue.jpg')); 
    cd(PATH_TO_ALL_USEFUL_PACKAGE);

    % Computation of the uncorrected p-value:
    t = tpdf(abs(slm_EffectsOnGroup.t), slm_EffectsOnGroup.df);
    uncorrected_pValues = double(t<=0.05).*t+double(t>0.05);
    clearvars struct; struct.P = uncorrected_pValues; struct.mask = mask; struct.thresh = thresholdUncorrectedPValue;
    SurfStatView( struct, averageSurface, [ '(ContrastNe-Uncorrected P-values )' Label1 '-' Label2]);
    cd(result_dir); 
    save2jpeg(strcat('Figures/', fileName, '/ContrastNegative-UncorrectedPValue.jpg')); 
    cd(PATH_TO_ALL_USEFUL_PACKAGE);

    % Computation of the corrected p-values: P-value threshold or statistic threshold for defining clusters, 0.001 by default
    [ pval, peak, clus ] = SurfStatP( slm_EffectsOnGroup , mask, clusterThreshold);
    pval.thresh = thresholdCorrectedPValue;
    SurfStatView( pval, averageSurface, ['(ContrastNe-Corrected P-values )' Label1 '-' Label2 ' (clusterThreshold = ' num2str(clusterThreshold) ')']);
    cd(result_dir); 
    save2jpeg(strcat('Figures/', fileName, '/ContrastNegative-CorrectedPValue.jpg')); 
    cd(PATH_TO_ALL_USEFUL_PACKAGE);
    
    disp('###')
    disp('After correction(Clusterwise Correction for Multiple Comparisons): ')
    if isempty(clus) ~= 1
        disp(['#Clusters found:                             ' num2str(  length(clus.P)        ) ])
        disp(['#Significative clusters (after correction) : ' num2str(  length(find(clus.P<=thresholdCorrectedPValue))  ) ])
    else
        disp(['No cluster found!'])
    end
    disp('###')
    
    % Computation of the false discovery rate :
    qval = SurfStatQ( slm_EffectsOnGroup , mask );
    SurfStatView( qval, averageSurface, ['ContrastNe-False discovery rate ' Label1 '-' Label2 ]);
    cd(result_dir);
    save2jpeg(strcat('Figures/', fileName, '/ContrastNegative-FalseDiscoveryRate.jpg' )); 
    cd(PATH_TO_ALL_USEFUL_PACKAGE);

