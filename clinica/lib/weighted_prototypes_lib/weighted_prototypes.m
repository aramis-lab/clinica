% Compute Streamline Weighted Prototypes of a bundle in VTK format
%
% Usage: weighted_prototypes(working_dir,filename_bundle,lambda_g,lambda_a,lambda_b,path_matlab_functions,path_CPP_code,path_Community_latest,bound_limit_input,degree_precision_input,num_iter_modularity_input,minimum_number_fibers_cluster_input,minValueTau_input,increase_radius_input)
%
% MANDATORY INPUTS:
% - working_dir: directory where files are saved
% - filename_bundle: filename of the fiber bundle which must be a .vtk file.
% It must have a list of 3D points and then it should use te keyword
% "LINES" for polygons. Each row describes a streamline. The first number
% is the number of points. The other numbers are the indexes of the points
% previously listed.
% - lambda_g: geometric kernel bandidth (as in usual currents)
% - lambda_a: kernel bandwidth of the STARTING structure
% - lambda_b: kernel bandwidth of the ENDING structure
% - path_matlab_functions
% - path_cpp_code: absolute path to C++ code 
% - path_community_latest: absolute path to Community detection folder
% To note: streamlines must have a consistent orientation, namely they must
% have the same starting and ending ROIs
% 
% OPTIONAL INPUTS:
% - bound_limit_input: maximum average angle (in radians) that a streamline 
% may have with the other streamlines in the framework of weighted currents. 
% Default value is 1.5359 = 88 degrees
% - degree_precision_input: percentage of the norm of the bundle explained by the
% weighted prototypes. Default value is 0.15, which means that the weighted
% prototypes will explain (1-0.15)*100 % of the norm of the bundle in the
% framework of weighted currents.
% - num_iter_modularity_input: Modularity computation is based on a greedy
% approach. Results may differ between iterations. The greater number of 
% iterations, the better. Default value is 10
% See "Fast unfolding of community hier archies in large networks", V. Blondel et al.
% - minimum_number_fibers_cluster_input: Clustering based on modularity may
% result in unbalanced clusters. We remove the clusters which have less
% than minimum_number_fibers_cluster_input fibers. Default value is 10
% - minValueTau_input: We remove the prototypes that approximate less than
% minValueTau_input fibers. Default value is 1
% - increase_radius_input: All tubes are normalised such that the maximum 
% radius is equal to 1mm. We then augment all radii of
% increase_radius_input. Default value is 0.02
%
% This function requires:
% - The binary files of the C++ functions in the folder cpp_code
% - CMake > 2.8
% - VTK > 6
% - ITK
% - Louvain community detection (https://sites.google.com/site/findcommunities/newversion/community.tgz?attredirects=0)
% which is already present 
% - Eigen (http://eigen.tuxfamily.org/index.php?title=Main_Page)
% which is also already present
%
%  Copyright Pietro GORI, Inria 
%  Written 16/08/2016
    
function [] = weighted_prototypes(working_dir,filename_bundle,lambda_g,lambda_a,lambda_b,path_matlab_functions,path_cpp_code,path_community_latest,bound_limit_input,degree_precision_input,num_iter_modularity_input,minimum_number_fibers_cluster_input,minValueTau_input,increase_radius_input)
 
    %% Check input parameters    
    switch nargin
        
        case 8            
            bound_limit=1.5359; % 85=1.4835 , 86=1.5010 , 87=1.5184 , 87.5=1.5272 , 88=1.5359 , 89=1.5533
            degree_precision=0.15; % it will explain (1-degree_precision)*100 % of the norm of the bundle                     
            num_iter_modularity=10;  
            minimum_number_fibers_cluster=10;
            minValueTau=1; 
            increase_radius=0.02; 
            
        case 14  
                         
            if bound_limit_input==0
                bound_limit=1.5359; 
            else
                bound_limit=bound_limit_input;
            end      
    
            if degree_precision_input==0
                degree_precision=0.15; 
            else
                degree_precision=degree_precision_input;
            end
            
            if num_iter_modularity_input==0
                num_iter_modularity=10;  
            else
                num_iter_modularity=num_iter_modularity_input;
            end
            
            if minimum_number_fibers_cluster_input==0
                minimum_number_fibers_cluster=10;
            else
                minimum_number_fibers_cluster=minimum_number_fibers_cluster_input;
            end

            if minValueTau_input==0
                minValueTau=1;
            else
                minValueTau=minValueTau_input;
            end            
   
            if increase_radius_input==0
                increase_radius=0.02;
            else
                increase_radius=increase_radius_input;
            end            
     
        otherwise
            error('Please insert either only working_dir,filename_bundle,lambda_g,lambda_a,lambda_b,path_matlab_functions,path_cpp_code,path_community_latest or all parameters. If you want to use default parameter please use 0 as input value')
            
    end   
    
    %% Matlab dependencies
    addpath(path_matlab_functions)
   
    %% Working directory
    cd(working_dir)
   
    %% Check dependencies
    if ~exist([ path_cpp_code '/bin/gramiam'],'file') 
        error(['Compile C++ code in the folder ' path_cpp_code '/bin inside cpp_code. Go to folder bin, type cmake ../ and then make. ITK and VTK are mandatory.'])
    end
    if ~exist([ path_cpp_code '/bin/medoids'],'file') 
        error(['Compile C++ code in the folder ' path_cpp_code '/bin inside cpp_code. Go to folder bin, type cmake ../ and then make. ITK and VTK are mandatory.'])
    end
    if ~exist([ path_cpp_code '/bin/write_tube'],'file') 
        error(['Compile C++ code in the folder ' path_cpp_code '/bin inside cpp_code. Go to folder bin, type cmake ../ and then make. ITK and VTK are mandatory.'])
    end
    if ~exist([ path_community_latest '/community'],'file') 
        error(['Compile C++ code in the folder community_latest. Just type make in the folder ' path_community_latest ])
    end
    if ~exist([ path_community_latest '/hierarchy'],'file') 
        error(['Compile C++ code in the folder community_latest. Just type make in the folder ' path_community_latest ])
    end

    disp(' ')
    disp(['Analysing bundle ' filename_bundle ' with parameters lambda_g: ' num2str(lambda_g) ' lambda_a:' num2str(lambda_a) ' lambda_b' num2str(lambda_b) ])
    disp(['Parameters: num_iter_modularity: ' num2str(num_iter_modularity) ', minValueTau: ' num2str(minValueTau) ', degree_precision: ' num2str(degree_precision) ', bound_limit: ' num2str(bound_limit) ])
    
    %% Parameters
%     filename_modularity='Modularity.mat';
    filename_vtk_no_outlier='NoOutliers.vtk';
    filename_vtk_clusters='Clusters.vtk';
    filename_medoids_tubes='Prototypes_tubes.vtk';
    filename_medoids_polyline='Prototypes.vtk';
    
    %% Loading bundle
    [Points,number_points_curve,] = vtk_bundles_polyline_reader(filename_bundle);

    %% Gramiam
    diary('gramiam.log')
    diary on                   
    eval(sprintf(['!' path_cpp_code '/bin/gramiam %s %i %f %f %f'],filename_bundle, 3, lambda_g, lambda_a, lambda_b))
    pause(1);
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
        eval(sprintf(['! ' path_community_latest '/community ' filename_graph_bin ' -l -1 -v -w ' filename_graph_weights ' > graph.tree']));
        diary off     

        fid = fopen('diary', 'r');
        if(fid==-1)
            error('Error: file descriptor not valid, check the file name');    
        end
        [~,~] = go_to_key_word(fid, 'Total');    
        line = fgetl(fid);
        Q = sscanf(line,'%f');
        fclose(fid);
        %disp(['Q: ' num2str(Q) ])
        delete('diary')     

        pause(1);

        eval(sprintf(['! ' path_community_latest '/hierarchy graph.tree > hierarchy.log']));    

        fid = fopen('hierarchy.log', 'r');
        if(fid==-1)
            error('Error: file descriptor not valid, check the file name');    
        end
        tline = fgetl(fid);
        num_level = sscanf(tline,'%*s %*s %*s %u');
        fclose(fid);

        %disp(['Number of levels: ' num2str(num_level) ])

        eval(sprintf(['! ' path_community_latest '/hierarchy graph.tree -l %u > C.log'], num_level-1));

        pause(1);

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
       Radius(i)=length(Fibers); % we consider the prototype
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

%     disp('Saving Modularity, NoOutliers and Clusters')
%     save(filename_modularity, 'Ci', 'Q_max', 'Cluster_no_outlier', '-v7.3');

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
    eval(sprintf(['! ' path_cpp_code '/bin/medoids %s %s %s %f %f %f %d %s'], filename_graph_bin, filename_graph_weights, filename_graph_diag, minValueTau, degree_precision, bound_limit, NClusters, filename_fibers_tot));
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

    %% Write tubes
    % Max radius is 1 otherwise problem to visualize

    Radius=TauGlobal;
    Radius=Radius/max(Radius);
    Radius=Radius+increase_radius; % if tubes are too small, increase this value         

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

    eval(sprintf(['! ' path_cpp_code '/bin/write_tube %s %s %s %s'], 'Points.txt','Number_Points_Curve_Medoids.txt', 'Radius.txt', filename_medoids_tubes));
        
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

    write_vtk_bundles_polyline(Points_Medoids,Number_points_curve_Medoids,Scalars_Medoids, [], [],filename_medoids_polyline)  

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

    write_vtk_bundles_polyline(Points_finale,number_points_curve_final,[],[],[],filename_vtk_no_outlier)
    write_vtk_bundles_polyline(Points_finale,number_points_curve_final,Scalars_color,[],[],filename_vtk_clusters)
    
    delete('C.log', 'graph.tree', 'index_fibers*', 'hierarchy.log', 'Medoids_global_normalised', 'Number_Points_Curve_Medoids.txt', 'Outliers_global', 'Points.txt', 'Radius.txt', 'Tau_global_normalised')
end
                    
                 


                