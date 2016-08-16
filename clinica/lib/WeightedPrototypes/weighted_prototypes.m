% Compute Streamline Weighted Prototypes of a bundle in VTK format
%
% Usage: weighted_prototypes(filename_bundle,lambda_g,lambda_a,lambda_b,bound_limit_input,degree_precision_input,num_iter_modularity_input,minimum_number_fibers_cluster_input,minValueTau_input)
%
% MANDATORY INPUTS:
% - filename_bundle:  the output file from recon-all pipeline,specifically, files: ?h.thickness.fwhm**.mgh.
% - lambda_g: the directory to contain the result images.
% - lambda_a: string, the linear model that fit into the GLM, for example '1+Lable'
% - lambda_b: string, the contrast that you want to use in the GLM, but the contrast should be inclued into the linearmodel, otherwise, you will get errors.
% 
% OPTIONAL INPUTS:
% - bound_limit_input: 
% - degree_precision_input: 
% - num_iter_modularity_input: 
% - minimum_number_fibers_cluster_input: 
% - minValueTau_input: 
%
%  Copyright Pietro GORI, Inria 
%  Written 16/08/2016
    
function [] = weighted_prototypes(filename_bundle,lambda_g,lambda_a,lambda_b,bound_limit_input,degree_precision_input,num_iter_modularity_input,minimum_number_fibers_cluster_input,minValueTau_input)
    
    addpath('Matlab_Functions')
    
    switch nargin
        
        case 4            
            bound_limit=1.5533; % 85=1.4835 , 86=1.5010 , 87=1.5184 , 87.5=1.5272 , 88=1.5359 , 89=1.5533
            degree_precision=0.15; % it will explain (1-degree_precision)*100 % of the total norm                     
            num_iter_modularity=5;  
            minimum_number_fibers_cluster=10;
            minValueTau=1;      
            
        case 9  
            
            if isempty(bound_limit_input)
                bound_limit=1.5533; % 85=1.4835 , 86=1.5010 , 87=1.5184 , 87.5=1.5272 , 88=1.5359 , 89=1.5533
            else
                bound_limit=bound_limit_input;
            end
            
            if isempty(degree_precision_input)
                degree_precision=0.15; % it will explain (1-degree_precision)*100 % of the total norm 
            else
                degree_precision=degree_precision_input;
            end
            
            if isempty(num_iter_modularity_input)
                num_iter_modularity=5;  
            else
                num_iter_modularity=num_iter_modularity_input;
            end
            
            if isempty(minimum_number_fibers_cluster_input)
                minimum_number_fibers_cluster=10;
            else
                minimum_number_fibers_cluster=minimum_number_fibers_cluster_input;
            end
            
            if isempty(minValueTau_input)
                minValueTau=1;
            else
                minValueTau=minValueTau_input;
            end 
     
        otherwise
            error('Please insert either only lambda_g,lambda_a,lambda_b and filename_bundle or all parameters. If you want to use default parameter please use [] as input value')
            
    end        

    disp(['Analysing bundle ' filename_bundle ' with parameters \lambda_g: ' num2str(lambda_g) ' \lambda_a:' num2str(lambda_a) ' \lambda_b' num2str(lambda_b) ])
    disp(['Parameters: num_iter_modularity: ' num2str(num_iter_modularity) ', minValueTau: ' num2str(minValueTau) ', degree_precision: ' num2str(degree_precision) ', bound_limit: ' num2str(bound_limit) ])
    
    filename_modularity=['Modularity_' filename_bundle(1:end-4) '.mat'];
    filename_vtk_NoOutlier=['NoOutliers_' filename_bundle(1:end-4) '.vtk'];
    filename_vtk_clusters=['Clusters_' filename_bundle(1:end-4) '.vtk'];
    filename_medoids_tubes=['Medoids_Tubes_' filename_bundle(1:end-4) '.vtk'];
    filename_Medoids_polyline=['Medoids_polyline_' filename_bundle(1:end-4) '.vtk'];
    
    %% Loading bundle
    [Points,number_points_curve,] = VTK_Bundles_Polyline_Reader(filename_bundle);

    %% Gramiam
    diary('Gramiam.log')
    diary on                   
    eval(sprintf('! CPP_code/Gramiam %s %i %f %f %f',filename_bundle, 3, lambda_g, lambda_a, lambda_b))
    pause(2);
    diary off

    filename_graph_bin='graph.bin';
    filename_graph_weights='graph.weights';
    filename_graph_diag='graph.diag';
    
    %% Modularity
    Data_C=cell(num_iter_modularity,1);
    Data_Q=zeros(num_iter_modularity,1);

    for i=1:num_iter_modularity

        diary('diary')
        diary on
        disp('Community construction')   
        eval(sprintf(['! Community_latest/community ' filename_graph_bin ' -l -1 -v -w ' filename_graph_weights ' > graph.tree']));
        diary off     

        fid = fopen('diary', 'r');
        if(fid==-1)
            error('Error: file descriptor not valid, check the file name');    
        end
        [~,~] = GoToKeyWord(fid, 'Total');    
        line = fgetl(fid);
        Q = sscanf(line,'%f');
        fclose(fid);
        %disp(['Q: ' num2str(Q) ])
        delete('diary')     

        pause(2);

        eval(sprintf('! Community_latest/hierarchy graph.tree > hierarchy.log'));    

        fid = fopen('hierarchy.log', 'r');
        if(fid==-1)
            error('Error: file descriptor not valid, check the file name');    
        end
        tline = fgetl(fid);
        num_level = sscanf(tline,'%*s %*s %*s %u');
        fclose(fid);

        %disp(['Number of levels: ' num2str(num_level) ])

        eval(sprintf('! Community_latest/hierarchy graph.tree -l %u > C.log', num_level-1));

        pause(2);

        C = load('C.log', '-ascii'); % C starts from 0 !!!
        C=C+1;  

        Data_C{i}=C(:,2);
        Data_Q(i)=Q;
    end

    [~,I_max]=max(Data_Q);
    Ci=Data_C{I_max};
    Number_clusters=length(unique(Ci));
    Groups=unique(Ci);

    disp('')
    disp(['Number clusters BEFORE number fiber constraint: ' num2str(Number_clusters) ])

    %% Creating Clusters
    Cluster=struct('Medoid',[],'Fibers',[],'Radius',[]);
    Radius=zeros(Number_clusters,1);
    index_to_remove=[];
    fibers_to_remove=[];

    for i=1:Number_clusters
       [Fibers,~]=find(Ci==Groups(i)); 
       Cluster(i).Fibers=Fibers;
       Radius(i)=length(Fibers); % I do consider the Medoid
       Cluster(i).Radius=Radius;
       if Radius(i)<minimum_number_fibers_cluster    
           index_to_remove=[index_to_remove ; i];  
           fibers_to_remove=[fibers_to_remove; Fibers];
       end
    end  

    %% Removing Outliers
    Cluster_no_outlier=Cluster;
    clear Cluster;
    Cluster_no_outlier(index_to_remove)=[];
    Radius(index_to_remove)=[];

    disp(['Minimum number fibers kept: ' num2str(minimum_number_fibers_cluster) ])
    disp(['Number clusters AFTER number fiber constraint: ' num2str(length(Radius)) ])
    disp('')

    disp('Saving Modularity, NoOutliers and Clusters')
    save(filename_modularity, 'Ci', 'Q_max', 'Cluster_no_outlier', '-v7.3');

    %% Prototypes  
      
    filename_final_med='Medoids_global_normalised';
    filename_final_tau='Tau_global_normalised';
    filename_final_out='Outliers_global';

    for i=1:length(Cluster_no_outlier)

        index_fibers=Cluster_no_outlier(i).Fibers;     

        filename_fibers=['index_fibers' num2str(i) ];
        fid = fopen(filename_fibers, 'w');                
        fwrite(fid,length(index_fibers),'int');  
        fwrite(fid,index_fibers-1,'int'); % In Matlab it starts from 1, in C++ from 0   
        fclose(fid);
    end

    NClusters = length(Cluster_no_outlier);
    filename_fibers_tot = 'index_fibers*';

    diary('Medoids.log')
    diary on
    disp('Prototypes')            
    eval(sprintf('! CPP_code/MedoidsFinale %s %s %s %f %f %f %d %s', filename_graph_bin, filename_graph_weights, filename_graph_diag, minValueTau, degree_precision, bound_limit, NClusters, filename_fibers_tot));
    pause(1);   
    diary off

            
    % SAVING RESULTS   
    fid = fopen(filename_final_med, 'r');
    MedoidsGlobal = fread(fid,'int');      
    MedoidsGlobal=MedoidsGlobal+1;
    fclose(fid);

    fid = fopen(filename_final_out, 'r');
    OutliersGlobal = fread(fid,'int');      
    OutliersGlobal=OutliersGlobal+1;
    fclose(fid);

    fid = fopen(filename_final_tau, 'r');
    TauGlobal = fread(fid,'single');  
    fclose(fid);

    Points_Medoids=[];
    Number_points_curve_Medoids=[];
    Points_Basal_Medoids=[];
    Points_Cortical_Medoids=[];

    for i=1:length(MedoidsGlobal)

        index=MedoidsGlobal(i); 
        Points_fibers=Points(sum(number_points_curve(1:(index-1)))+1:sum(number_points_curve(1:index)),:);

        if size(Points_fibers,1)~=number_points_curve(index)
            error('PROBLEM! NUMBER POINTS DO NOT CORRESPOND')
        end  

        Points_Medoids=[Points_Medoids; Points_fibers];  
        Points_Basal_Medoids=[Points_Basal_Medoids ; Points_fibers(1,:)];
        Points_Cortical_Medoids=[Points_Cortical_Medoids ; Points_fibers(end,:)];
        Number_points_curve_Medoids=[Number_points_curve_Medoids ; number_points_curve(index)];

    end
    
    disp('Writing Tube')

    %% Max radius is 1 otherwise problem to visualize

    Radius=TauGlobal;
    Radius=Radius/max(Radius);
    Radius=Radius+0.2;         

    % Points Medoids
    filename_points = 'Points.txt';
    fid = fopen(filename_points, 'w'); 
    nPts=size(Points_Medoids,1);
    for i=1:nPts
        fprintf(fid, '%f %f %f\n',Points_Medoids(i,1),Points_Medoids(i,2),Points_Medoids(i,3));
    end
    fclose(fid);

    % Number_points_curve_Medoids
    filename_Number_Points_Curve_Medoids = 'Number_Points_Curve_Medoids.txt';
    fid = fopen(filename_Number_Points_Curve_Medoids, 'w'); 
    nMedoid=size(Number_points_curve_Medoids,1);
    for i=1:nMedoid
        fprintf(fid, '%d\n',Number_points_curve_Medoids(i));
    end	
    fclose(fid);

    % Radius
    filename_Radius = 'Radius.txt';
    fid = fopen(filename_Radius, 'w'); 
    if nMedoid~=size(Radius,1)
        disp('PROBLEM')
    end
    for i=1:nMedoid
        fprintf(fid, '%d\n',Radius(i));
    end		
    fclose(fid);

    eval(sprintf('! CPP_code/WriteTube %s %s %s %s', 'Points.txt','Number_Points_Curve_Medoids.txt', 'Radius.txt', filename_medoids_tubes));
        
    % Writing Prototypes    
    Scalars_Medoids=zeros(size(Points_Medoids,1),1);
    count=1;
    for i=1:length(Number_points_curve_Medoids)    
        for j=1:Number_points_curve_Medoids(i)        
            Scalars_Medoids(count)=TauGlobal(i);
            count=count+1;
        end    
    end

    if count-1 ~= size(Points_Medoids,1)
        error('Error')
    end

    Write_vtk_bundles_polyline(Points_Medoids,Number_points_curve_Medoids,Scalars_Medoids, [], [],filename_Medoids_polyline)  

    % Writing NoOutliers and Clusters   
    Points_Basal=[];
    Points_Cortical=[];
    Index_Points=[];

    numout=0;
    number_points_curve_final=[];
    Points_finale=[];

    Scalars_color=[];
    colors=linspace(1,length(Cluster_no_outlier),length(Cluster_no_outlier));
    idx=randperm(length(Cluster_no_outlier));
    colors=colors(idx);

    for i=1:length(Cluster_no_outlier)
        index_fibers=Cluster_no_outlier(i).Fibers; 
        for j=1:length(index_fibers) 
            index=index_fibers(j);
            test=intersect(OutliersGlobal,index);
            if isempty(test)
                Points_fibers=Points(sum(number_points_curve(1:(index-1)))+1:sum(number_points_curve(1:index)),:);
                Points_finale=[Points_finale; Points_fibers];
                number_points_curve_final=[number_points_curve_final; size(Points_fibers,1) ];
                Scalars_color=[Scalars_color; ones(size(Points_fibers,1),1)*colors(i)];
                Points_Basal=[Points_Basal ; Points_fibers(1,:)];
                Points_Cortical=[Points_Cortical ; Points_fibers(end,:)];
                Index_Points=[Index_Points; index];
            else
                numout=numout+1;
            end
        end  
    end

    if (numout ~= length(OutliersGlobal))
        error('Problem with outliers!')
    end  

    Write_vtk_bundles_polyline(Points_finale,number_points_curve_final,[],[],[],filename_vtk_NoOutlier)
    Write_vtk_bundles_polyline(Points_finale,number_points_curve_final,Scalars_color,[],[],filename_vtk_clusters)
    
end
                    
                 


                