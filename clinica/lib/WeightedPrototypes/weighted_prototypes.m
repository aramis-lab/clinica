function [] = weighted_prototypes(lambda_g,lambda_a,lambda_b,filename_bundle)
                    
    bound_limit=1.5533; % 85=1.4835 , 86=1.5010 , 87=1.5184 , 87.5=1.5272 , 88=1.5359 , 89=1.5533
    degree_precision=0.15; % it will explain (1-degree_precision)*100 % of the total norm                     
    num_iter_modularity=5;  
    minimum_number_fibers_cluster=10;
    minValueTau=1;
    dim=3;

%     max_n_fibers=100000;

    disp(['Analysing bundle ' filename_bundle ' with parameters \lambda_g: ' num2str(lambda_g) ' \lambda_a:' num2str(lambda_a) ' \lambda_b' num2str(lambda_b) ])
    filename_gram='graph.bin';
    filename_modularity=['Modularity_' kind '_' type '_full_W' num2str(W) '_C' num2str(Fc) '_B' num2str(Fb) '.mat'];
    filename_vtk_NoOutlier=['NoOutliers_' kind '_' type '_full_W' num2str(W) '_C' num2str(Fc) '_B' num2str(Fb) '.vtk'];
    filename_vtk_clusters=['Clusters_' kind '_' type '_full_W' num2str(W) '_C' num2str(Fc) '_B' num2str(Fb) '.vtk'];
    filename_medoids_tubes=['Medoids_Tubes_' kind '_' type '_full_W' num2str(W) '_C' num2str(Fc) '_B' num2str(Fb) '.vtk'];
    filename_Medoids_polyline=['Medoids_polyline_' kind '_' type '_full_W' num2str(W) '_C' num2str(Fc) '_B' num2str(Fb) '.vtk'];

%     filename_result=['Result_' kind '_' type '_full_W' num2str(W) '_C' num2str(Fc) '_B' num2str(Fb) '.mat' ];
%     density_fibers_basal=['Density_fibers_basal_' kind '_' type '_full_W' num2str(W) '_C' num2str(Fc) '_B' num2str(Fb) '.vtk'];
%     density_med_basal=['Density_med_basal_' kind '_' type '_full_W' num2str(W) '_C' num2str(Fc) '_B' num2str(Fb) '.vtk'];
%     density_fibers_cortex=['Density_fibers_cortex_' kind '_' type '_full_W' num2str(W) '_C' num2str(Fc) '_B' num2str(Fb) '.vtk'];
%     density_med_cortex=['Density_med_cortex_' kind '_' type '_full_W' num2str(W) '_C' num2str(Fc) '_B' num2str(Fb) '.vtk'];
%     filename_density_uniform_B=['Density_uniform_basal_' kind '_' type '_full_W' num2str(W) '_C' num2str(Fc) '_B' num2str(Fb) '.vtk'];
%     filename_density_uniform_C=['Density_uniform_cortex_' kind '_' type '_full_W' num2str(W) '_C' num2str(Fc) '_B' num2str(Fb) '.vtk'];
%     filename_fibers_uniform_down=['Uniform_' kind '_' type '_full_W' num2str(W) '_C' num2str(Fc) '_B' num2str(Fb) '.vtk'];



    disp(['Parameters: num_iter_modularity: ' num2str(num_iter_modularity) ', minValueTau: ' num2str(minValueTau) ', degree_precision: ' num2str(degree_precision) ', bound_limit: ' num2str(bound_limit) ])


    %% Downsampling
%     [Points,number_points_curve] = VTK_Bundles_Polyline_Reader(filename_bundle);          
%     Number_fibers_original=length(number_points_curve);
%     if (Number_fibers_original>max_n_fibers)
%         if fopen([ filename_bundle(1:end-9) '_down_' num2str(max_n_fibers) '.vtk' ])==(-1)                 
%             filename=filename_bundle;
%             new_filename=[ filename_bundle(1:end-9) '_down.vtk' ];        
%             [New_Points,New_number_points_curve] = DownSampleBundlesVTK_uniform_number_fast(filename,new_filename,max_n_fibers,Points,number_points_curve); 
%             Points=New_Points;
%             number_points_curve=New_number_points_curve;
%             filename_bundle=[ filename_bundle(1:end-9) '_down_' num2str(max_n_fibers) '.vtk' ];
%         else                    
%             [Points,number_points_curve] = VTK_Bundles_Polyline_Reader([ filename_bundle(1:end-9) '_down_' num2str(max_n_fibers) '.vtk' ]);
%             filename_bundle=[ filename_bundle(1:end-9) '_down_poly_' num2str(max_n_fibers) '.vtk' ];
%         end
% 
%         Number_fibers_original=length(number_points_curve);
% 
%         if (Number_fibers_original ~= max_n_fibers)
%             error('PROBLEM WITH DOWNSAMPLING')
%         end 
%     end         



          

            if fopen(filename_gram)==(-1)

                if ~isempty(list)

                    if (length(list)==1)
                        if (strcmp(filename_dossier,list(1).name))
                            error('DELETE folder and start again')
                        end
                    end      

                    disp('Using Gram matrix already computed')
                    disp(' ')
                    trovato=0;                
                    for tt=1:length(list)
                        if (~strcmp(filename_dossier,list(tt).name))
                            if fopen([list(tt).name '/' filename_gram])~=(-1)
                                filename_graph_bin=['../' list(tt).name '/graph.bin'];
                                filename_graph_weights=['../' list(tt).name '/graph.weights'];
                                filename_graph_diag=['../' list(tt).name '/graph.diag'];
                                trovato=1;
                                list_mod=dir([list(tt).name '/Modularity' '*']);
                                if ~isempty(list_mod) && length(list_mod)==1
                                    load([list(tt).name '/' list_mod.name ])
                                    trovato_mod=1;
                                end
                                break; 
                            end
                        end            
                    end

                    if trovato==0
                        error('There is a problem.. delete all folders and start again')
                    end

                    cd(filename_dossier)  

                else

                    cd(filename_dossier)  
                    diary('Gramiam.log')
                    diary on
                    if strcmp(type,'det')
                        eval(sprintf(['! ~/Cpp/Prototypes/' computer '/Gramiam %s %i %f %f %f'],['../../Fiber_bundles_det_1s/' filename_bundle], dim, W, Fb, Fc))
                        %eval(sprintf(['! ~/Cpp/Gram_matrix_Fiber_Bundle/' computer '/GramiamOLD %s %i %f %f %f %s %s'],['../../Fiber_bundles_det_1s/' filename_bundle], dim, W, Fb, Fc, filename_gram, filename_diag))
                    elseif strcmp(type,'prob')
                        eval(sprintf(['! ~/Cpp/Prototypes/' computer '/Gramiam %s %i %f %f %f'],['../../Fiber_bundles_prob_8s/' filename_bundle], dim, W, Fb, Fc))
                        %eval(sprintf(['! ~/Cpp/Gram_matrix_Fiber_Bundle/' computer '/GramiamOLD %s %i %f %f %f %s %s'],['../../Fiber_bundles_prob_8s/' filename_bundle], dim, W, Fb, Fc, filename_gram, filename_diag))
                    else
                        error('ERROR! DET OR PROB!')
                    end
                    pause(2);
                    diary off
                    
                    filename_graph_bin='graph.bin';
                    filename_graph_weights='graph.weights';
                    filename_graph_diag='graph.diag';

                end

            else

                cd(filename_dossier)  
                filename_graph_bin='graph.bin';
                filename_graph_weights='graph.weights';
                filename_graph_diag='graph.diag';
                disp('Gramiam Matrix already computed')
                disp(' ')
            end


            %% Modularity           
            if fopen(filename_modularity)==(-1) && (trovato_mod==0)

                %disp('Converting to binary file')
                %eval(sprintf(['! /lena13/home_users/users/gori/Cpp/Community_latest_' computer '/convert -i %s -o graph.bin -w graph.weights'], filename_gram))
                %pause(2);

                Data_C=cell(num_iter_modularity,1);
                Data_Q=zeros(num_iter_modularity,1);

                for i=1:num_iter_modularity

                    diary('diary')
                    diary on
                    disp('Community construction')   
                    eval(sprintf(['! ~/Cpp/Community_latest/' computer '/community ' filename_graph_bin ' -l -1 -v -w ' filename_graph_weights ' > graph.tree']));
                    diary off     

                    fid = fopen('diary', 'r');
                    if(fid==-1)
                        error('Error: file descriptor not valid, check the file name');    
                    end
                    [out,LINE] = GoToKeyWord(fid, 'Total');    
                    line = fgetl(fid);
                    Q = sscanf(line,'%f');
                    fclose(fid);
                    %disp(['Q: ' num2str(Q) ])
                    delete('diary')     

                    pause(2);

                    eval(sprintf(['! ~/Cpp/Community_latest/' computer '/hierarchy graph.tree > hierarchy.log']));    

                    fid = fopen('hierarchy.log', 'r');
                    if(fid==-1)
                        error('Error: file descriptor not valid, check the file name');    
                    end
                    tline = fgetl(fid);
                    num_level = sscanf(tline,'%*s %*s %*s %u');
                    fclose(fid);

                    disp(['Number of levels: ' num2str(num_level) ])

                    eval(sprintf(['! ~/Cpp/Community_latest/' computer '/hierarchy graph.tree -l %u > C.log'], num_level-1));

                    pause(2);

                    C = load('C.log', '-ascii'); % C starts from 0 !!!
                    C=C+1;  

                    Data_C{i}=C(:,2);
                    Data_Q(i)=Q;
                end

                %% ANALYSIS RESULTS

                [Q_max,I_max]=max(Data_Q);
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
                %        disp(['PROBLEM! CLUSTER ' num2str(i)  ' HAS LESS THAN ' num2str(minimum_number_fibers_cluster) ' FIBERS'])
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

    %             number_points_curve_no_outliers_color=number_points_curve;
    %             number_points_curve_no_outliers_color(fibers_to_remove)=[];

                disp('Saving Modularity, NoOutliers and Clusters')
                save(filename_modularity, 'Ci', 'Q_max', 'Cluster_no_outlier', '-v7.3');

            else

                disp('Modularity already computed. Loading it!')
                
                if (trovato_mod==0)
                    load(filename_modularity)
                end

            end

            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Medoids %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            filename_final_med='Medoids_global_normalised';
            filename_final_tau='Tau_global_normalised';
            filename_final_out='Outliers_global';

            if fopen(filename_final_out)==(-1) && ~strcmpi(computer,'blade')
  
                filename_medoids='Medoids_global';
                filename_tau='Tau_global';
                filename_outliers='Outliers_global';

                for i=1:length(Cluster_no_outlier)

                    index_fibers=Cluster_no_outlier(i).Fibers;     

                    filename_fibers=['index_fibers' num2str(i) ];
                    fid = fopen(filename_fibers, 'w');
                    %fprintf(fid,'%u\n',index_fibers-1'); % In Matlab it starts from 1, in C++ from 0
                    fwrite(fid,length(index_fibers),'int');  
                    fwrite(fid,index_fibers-1,'int'); % In Matlab it starts from 1, in C++ from 0   
                    fclose(fid);

                end

                NClusters = length(Cluster_no_outlier);
                filename_fibers_tot = 'index_fibers*';
                
                diary('Medoids.log')
                diary on
                
                disp('Medoids')
                % tic
                eval(sprintf(['! ~/Cpp/Prototypes/' computer '/MedoidsFinale %s %s %s %f %f %f %d %s'], filename_graph_bin, filename_graph_weights, filename_graph_diag, minValueTau, degree_precision, bound_limit, NClusters, filename_fibers_tot));
                % toc
                pause(1);   
                
                diary off

            elseif strcmpi(computer,'blade')
                disp('Medoids can not be computed on blade')
                
            else
                disp('Medoids already computed. Loading it!')
            end
            
            %% SAVING RESULTS     
            if ~strcmpi(computer,'blade')

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
            end

            if fopen(filename_medoids_tubes)==(-1) && ~strcmpi(computer,'blade')
                
                disp('Writing Tube')

                %% Max radius is 1 otherwise problem to visualize
                %copyfile('~/matlab/Tubes_VTK/Write_Tube','./Write_Tube');

                % Radius=Radius/max(Radius);
                % Radius=log(Radius);
                Radius=TauGlobal;
                Radius=Radius/max(Radius);
                Radius=Radius+0.1;
                % cd('/lena13/home_users/users/gori/matlab/Tubes_VTK')

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

                eval(sprintf(['! ~/Cpp/Prototypes/' computer '/WriteTube %s %s %s %s'], 'Points.txt','Number_Points_Curve_Medoids.txt', 'Radius.txt', filename_medoids_tubes));

            end

            %% Scalars needs to be same dimensiona as Points

            if fopen(filename_result)==(-1) && ~strcmpi(computer,'blade')
                
                disp('Writing Results.mat')

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

                % It should be for WriteTubeColor but it does not work...
                % % Colors
                % filename_Colors = 'Colors.txt';
                % fid = fopen(filename_Colors, 'w'); 
                % for i=1:size(Points_Medoids,1)
                %     fprintf(fid, '%d\n',Scalars(i));
                % end		
                % fclose(fid);

                save(filename_result,'Points_Medoids','Number_points_curve_Medoids','Scalars_Medoids','-v7.3')
            end

            disp(' ')
            %% Saving Densities Total
            
            if ~strcmpi(computer,'blade')
            
                load(filename_result)

                if fopen([filename_vtk_NoOutlier(1:end-4) '.mat'])==(-1) || exist('Index_Points','var')==0

                    disp('Writing No Outliers')

                    if (Number_fibers_original~=length(number_points_curve))
                        error('Problem! modified number_points_curve... ')
                    end

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

                    save([filename_vtk_NoOutlier(1:end-4) '.mat'],'Points_finale','number_points_curve_final','-v7.3')
                    Write_vtk_bundles_segments(Points_finale,number_points_curve_final,[],[],[],filename_vtk_NoOutlier)
                    Write_vtk_bundles_segments(Points_finale,number_points_curve_final,Scalars_color,[],[],filename_vtk_clusters)

                    save(filename_result,'Points_Medoids','Number_points_curve_Medoids','Scalars_Medoids','Points_finale','number_points_curve_final','Scalars_color','Index_Points','Points_Cortical','Points_Basal','-v7.3')
                 else
                     disp('Bundle without outlier and Cluster already computed')                 
                     disp(' ')
                 end

                if ( fopen(density_fibers_cortex)==(-1) && density==1)     

                    diary('density.log')
                    diary on

                    cd('../../cortical_surface')
                    if fopen([cortex(1:end-4) '.mat'])==(-1)    
                        [Vertex_CS,Faces_CS] = VTKPolyDataReader(cortex);
                        save([cortex(1:end-4) '.mat'],'Vertex_CS','Faces_CS','-v7.3');    
                    else
                        load([cortex(1:end-4) '.mat'])
                    end

                    cd('../sub_cortical')
                    if fopen([basal(1:end-4) '.mat'])==(-1)    
                        [Vertex_BG,Faces_BG] = VTKPolyDataReader(basal);
                        save([basal(1:end-4) '.mat'],'Vertex_BG','Faces_BG','-v7.3');    
                    else
                        load([basal(1:end-4) '.mat'])
                    end

                    cd(['../Prototypes/' filename_dossier ])   

                    %% Basal
                    disp('Density Basal')
                    [c, ia, ib] = intersect(Index_Points,MedoidsGlobal);

                    if (sum(c-MedoidsGlobal)~=0)
                        error('ERROR! Theu should be equal...')
                    end

                    if (sum(Index_Points(ia)-MedoidsGlobal)~=0)
                        error('ERROR! Theu should be equal...')
                    end  

                    if (sum(sum(Points_Cortical(ia,:)-Points_Cortical_Medoids(ib,:)))~=0) 
                        error('It should be 0')
                    end    

                    if length(ia)~= length(MedoidsGlobal)
                        erro('They should be equal!')
                    end

                    %     [M,~]=size(Vertex_BG);
                    NPbasal=size(Points_Basal,1);
                    NMedoids=size(Points_Basal_Medoids,1);

                %     if (M>20000 && NPbasal>20000)       

                        Points_Basal_reshape=reshape(Points_Basal',NPbasal*3,1);
                        filename_fibers='PointsBasal';
                        fid = fopen(filename_fibers, 'w');
                        fwrite(fid,NPbasal,'int');  
                        fwrite(fid,Points_Basal_reshape,'single');
                        fclose(fid);

                        Points_Medoids_reshape=reshape(Points_Basal_Medoids',NMedoids*3,1);
                        filename_medoids_basal='PointsBasalMedoids';
                        fid = fopen(filename_medoids_basal, 'w');
                        fwrite(fid,NMedoids,'int');  
                        fwrite(fid,Points_Medoids_reshape,'single');
                        fclose(fid);


                %         filenameIndex='IndexPointsMedoids';
                %         fid = fopen(filenameIndex, 'w');  
                %         fwrite(fid,length(ia),'int');  
                %         fwrite(fid,ia-1,'int'); % In Matlab it starts from 1, in C++ from 0   
                %         fclose(fid);

                        eval(sprintf(['! ~/Cpp/Prototypes/' computer '/Density %s %s %s %s %f'], filename_fibers, filename_medoids_basal, filename_final_tau, ['../../sub_cortical/' basal], lambda_B));

                        fid = fopen('DensityFibers.bin', 'r');
                        DensityFibers = fread(fid,'single');    
                        fclose(fid);

                        fid = fopen('DensityMedoids.bin', 'r');
                        DensityMedoids = fread(fid,'single');    
                        fclose(fid);

                %     else       
                % 
                %         Gamma_B=zeros(M,NPbasal);
                % 
                %         for i=1:3
                %            Matrix_V=repmat(Vertex_BG(:,i),1,NPbasal);
                %            Matrix_F=repmat(Points_Basal(:,i)',M,1);
                %            Matrix_Diff=(Matrix_V-Matrix_F).^2;
                %            Gamma_B=Gamma_B+Matrix_Diff;    
                %         end
                %         Gamma_B=exp(-Gamma_B./(2*lambda_B^2));
                %         Gamma_B=Gamma_B.*((1/(2*pi)^(3/2))*(1/lambda_B^3));
                % 
                %         DensityFibers=sum(Gamma_B,2)/NPbasal; % density of all points of a fiber on each vertex
                % 
                %         %% Compute density Prototypes Basal Ganglia
                %         Gamma_Med_B=Gamma_B(:,ia);
                %         if size(TauGlobal,2)==1
                %             TauGlobal=TauGlobal';
                %         end
                %         Gamma_Med_B=Gamma_Med_B.*repmat(TauGlobal,M,1);
                %         DensityMedoids=sum(Gamma_Med_B,2)/sum(TauGlobal); % density of all medoids on each vertex 
                % 
                %         disp(['N is: ' num2str(NPbasal) ' and sum tau is: ' num2str(sum(TauGlobal)) ])
                % 
                %     end

                    % Compute Test

                    Emp_distr_F_basal=cumsum(DensityFibers)./(sum(DensityFibers)); % In this way integral density is 1
                    Emp_distr_Med=cumsum(DensityMedoids)./(sum(DensityMedoids));

                    [D,I]=max(abs(Emp_distr_F_basal-Emp_distr_Med));

                    % http://www.jstor.org/stable/2280095?seq=1#page_scan_tab_contents 
                    % level alpha of 0.05
                    % dalpha=1.95/sqrt(sum(TauGlobal));
                    dalpha=1.36*sqrt((NPbasal+length(MedoidsGlobal))/(NPbasal*length(MedoidsGlobal)));

                    disp('Testing Proto=Fibers')
                    if D>dalpha
                        disp(['Hypothesis rejected! With a level alpha of 0.001, percentage of diff = ' num2str((D-dalpha)*100/dalpha) '%' ])
                    else
                        disp('Hypothesis acceted!')
                    end   

                    %% PLOT
                    if show==1
                        figure(2)
                        set(figure(2),'Position',[ 1921    1    1920   1121])
                        hold on
                        plot(1:length(Emp_distr_F_basal),Emp_distr_F_basal,'r',1:length(Emp_distr_F_basal),Emp_distr_Med,'k');
                        plot(1:length(Emp_distr_F_basal),Emp_distr_F_basal+dalpha,'m');
                        plot(1:length(Emp_distr_F_basal),Emp_distr_F_basal-dalpha,'m');
                        if D>dalpha
                            plot([I I],[0,1],'g');
                        end
                        legend('fibers','proto','accepted boundaries')
                        xlabel('X')
                        ylabel('Cumulative Probability')
                    end

                    %% Coloring 

                    Scalars_F = DensityFibers;
                    Scalars_F=Scalars_F./sum(Scalars_F);

                    Scalars_Proto = DensityMedoids;
                    Scalars_Proto=Scalars_Proto./sum(Scalars_Proto);

                    VTKPolyDataWriter(Vertex_BG, Faces_BG, Scalars_F, [], [], density_fibers_basal);
                    VTKPolyDataWriter(Vertex_BG, Faces_BG, Scalars_Proto, [], [], density_med_basal);

                    %[Scalars_F,Scalars_Proto] = density_fiber_surface_total(Points_Basal,Index_Points,MedoidsGlobal,TauGlobal,Vertex_BG,Faces_BG,lambda_B,1,1,'prova_fibers.vtk','prova_med.vtk');

                    %% Cortex
                    disp(' ')
                    disp('Density Cortex')

                %     [M,~]=size(Vertex_CS);
                    NPcortex=size(Points_Cortical,1);

                    if NPcortex>30000        
                        INDEX=randperm(NPcortex,30000); 
                        Points_Cortical=Points_Cortical(INDEX,:);
                        NPcortex=30000;
                    end

                %     if (M*NPcortex>245772000)       

                        Points_Cortex_reshape=reshape(Points_Cortical',NPcortex*3,1);
                        filename_fibers='PointsCortex';
                        fid = fopen(filename_fibers, 'w');
                        fwrite(fid,NPcortex,'int');  
                        fwrite(fid,Points_Cortex_reshape,'single');
                        fclose(fid);    

                        Points_Medoids_reshape=reshape(Points_Cortical_Medoids',NMedoids*3,1);
                        filename_medoids_cortex='PointsCortexMedoids';
                        fid = fopen(filename_medoids_cortex, 'w');
                        fwrite(fid,NMedoids,'int');  
                        fwrite(fid,Points_Medoids_reshape,'single');
                        fclose(fid);

                %         NIndex=length(ia);
                %         if NIndex~= length(MedoidsGlobal)
                %             erro('They should be equal!')
                %         end
                %         filenameIndex='IndexPointsMedoids';
                %         fid = fopen(filenameIndex, 'w');  
                %         fwrite(fid,NIndex,'int');  
                %         fwrite(fid,ia-1,'int'); % In Matlab it starts from 1, in C++ from 0   
                %         fclose(fid);

                        eval(sprintf(['! ~/Cpp/Prototypes/' computer '/Density %s %s %s %s %f'], filename_fibers, filename_medoids_cortex, filename_final_tau, ['../../cortical_surface/' cortex], lambda_C));

                        fid = fopen('DensityFibers.bin', 'r');
                        DensityFibers = fread(fid,'single');    
                        fclose(fid);

                        fid = fopen('DensityMedoids.bin', 'r');
                        DensityMedoids = fread(fid,'single');    
                        fclose(fid);

                %     else       
                % 
                %         Gamma_B=zeros(M,NPcortex);
                % 
                %         for i=1:3
                %            Matrix_V=repmat(Vertex_CS(:,i),1,NPcortex);
                %            Matrix_F=repmat(Points_Cortical(:,i)',M,1);
                %            Matrix_Diff=(Matrix_V-Matrix_F).^2;
                %            Gamma_B=Gamma_B+Matrix_Diff;    
                %         end
                %         Gamma_B=exp(-Gamma_B./(2*lambda_C^2));
                %         Gamma_B=Gamma_B.*((1/(2*pi)^(3/2))*(1/lambda_C^3));
                % 
                %         DensityFibers=sum(Gamma_B,2)/NPcortex; % density of all points of a fiber on each vertex
                % 
                %         %% Compute density Prototypes Basal Ganglia
                %         Gamma_Med_B=Gamma_B(:,ia);
                %         if size(TauGlobal,2)==1
                %             TauGlobal=TauGlobal';
                %         end
                %         Gamma_Med_B=Gamma_Med_B.*repmat(TauGlobal,M,1);
                %         DensityMedoids=sum(Gamma_Med_B,2)/sum(TauGlobal); % density of all medoids on each vertex 
                % 
                %         disp(['N is: ' num2str(NPcortex) ' and sum tau is: ' num2str(sum(TauGlobal)) ])
                % 
                %     end

                %     %% Prova (to delete)
                %     Points_Cortical_medoids=Points_Cortical(ia,:);
                %     Scalars = density_prototypes(Points_Cortical_medoids,TauGlobal,Vertex_CS,Faces_CS,lambda_C,0,1,'prova_density_cortex.vtk');

                    % Compute Test

                    Emp_distr_F_cortex=cumsum(DensityFibers)./(sum(DensityFibers)); % In this way integral density is 1
                    Emp_distr_Med=cumsum(DensityMedoids)./(sum(DensityMedoids));

                    [D,I]=max(abs(Emp_distr_F_cortex-Emp_distr_Med));

                    % http://www.jstor.org/stable/2280095?seq=1#page_scan_tab_contents 
                    % level alpha of 0.001  
                    dalpha=1.36*sqrt((NPcortex+length(MedoidsGlobal))/(NPcortex*length(MedoidsGlobal)));

                    disp('Testing Proto=Fibers')
                    if D>dalpha
                        disp(['Hypothesis rejected! With a level alpha of 0.001, percentage of diff = ' num2str((D-dalpha)*100/dalpha) '%' ])
                    else
                        disp('Hypothesis acceted!')
                    end   

                    %% PLOT
                    if show==1
                        figure(3)
                        set(figure(3),'Position',[ 1921    1    1920   1121])
                        hold on
                        plot(1:length(Emp_distr_F_cortex),Emp_distr_F_cortex,'r',1:length(Emp_distr_F_cortex),Emp_distr_Med,'k');
                        plot(1:length(Emp_distr_F_cortex),Emp_distr_F_cortex+dalpha,'m');
                        plot(1:length(Emp_distr_F_cortex),Emp_distr_F_cortex-dalpha,'m');
                        if D>dalpha
                            plot([I I],[0,1],'g');
                        end
                        legend('fibers','proto','accepted boundaries')
                        xlabel('X')
                        ylabel('Cumulative Probability')
                    %     bar(Density_F_B)
                    %     figure
                    %     bar(Density_Med_B) 
                    end

                    %% Coloring 

                    Scalars_F = DensityFibers;
                    Scalars_F=Scalars_F./sum(Scalars_F);

                    Scalars_Proto = DensityMedoids;
                    Scalars_Proto=Scalars_Proto./sum(Scalars_Proto);

                    VTKPolyDataWriter(Vertex_CS, Faces_CS, Scalars_F, [], [], density_fibers_cortex);
                    VTKPolyDataWriter(Vertex_CS, Faces_CS, Scalars_Proto, [], [], density_med_cortex);

                    % %[~,Scalars_Med_C] = density_fiber_surface_total(Points_Cortical,Index_Points,MedoidsGlobal,TauGlobal,Vertex_CS,Faces_CS,lambda_C,show,save_vtk,density_fibers_cortex,density_med_cortex);


                    %% UNIFORM DOWNSAMPLING
                    total_fibers=length(number_points_curve_final);
                    n_fibers=length(MedoidsGlobal); 
                    INDEX=randperm(total_fibers,n_fibers); 

                    Points_unif_down_B = zeros(length(INDEX),3);
                    Points_unif_down_C = zeros(length(INDEX),3);

                    New_number_points_curve_uniform=number_points_curve_final(INDEX);
                    New_Points_uniform=[];

                    for j=1:length(INDEX)
                        index=INDEX(j);
                        Points_fibers=Points_finale(sum(number_points_curve_final(1:(index-1)))+1:sum(number_points_curve_final(1:index)),:);
                        New_Points_uniform=[New_Points_uniform; Points_fibers];
                        Points_unif_down_B(j,:) = Points_fibers(1,:);
                        Points_unif_down_C(j,:) = Points_fibers(end,:);              
                    end

                    [Scalars_F_B] = density_fiber_surface_fibers(Points_unif_down_B,Vertex_BG,Faces_BG,lambda_B,0,1,filename_density_uniform_B);
                    [Scalars_F_C] = density_fiber_surface_fibers(Points_unif_down_C,Vertex_CS,Faces_CS,lambda_C,0,1,filename_density_uniform_C);

                    %% TEST
                    disp(' ')
                    disp('Density Fibers Vs Uniform')

                    Emp_distr_Uniform_basal=cumsum(Scalars_F_B)./(sum(Scalars_F_B));
                    Emp_distr_Uniform_cortex=cumsum(Scalars_F_C)./(sum(Scalars_F_C));  

                    % Compute Test
                    [D,I]=max(abs(Emp_distr_F_basal-Emp_distr_Uniform_basal));
                    dalpha=1.36*sqrt((total_fibers+n_fibers)/(total_fibers*n_fibers));

                    disp('Testing Uniform=Fibers Basal')
                    if D>dalpha
                        disp(['Hypothesis rejected! With a level alpha of 0.001, percentage of diff = ' num2str((D-dalpha)*100/dalpha) '%' ])
                    else
                        disp('Hypothesis acceted!')
                    end   
                    disp(' ')

                    if show==1
                        figure(4)
                        set(figure(4),'Position',[ 1921    1    1920   1121])
                        hold on
                        plot(1:length(Emp_distr_F_basal),Emp_distr_F_basal,'r',1:length(Emp_distr_F_basal),Emp_distr_Uniform_basal,'k');
                        plot(1:length(Emp_distr_F_basal),Emp_distr_F_basal+dalpha,'m');
                        plot(1:length(Emp_distr_F_basal),Emp_distr_F_basal-dalpha,'m');
                        if D>dalpha
                            plot([I I],[0,1],'g');
                        end
                        legend('fibers','proto','accepted boundaries')
                        xlabel('X')
                        ylabel('Cumulative Probability')
                    %     bar(Density_F_B)
                    %     figure
                    %     bar(Density_Med_B) 
                    end

                    % Compute Test
                    [D,I]=max(abs(Emp_distr_F_cortex-Emp_distr_Uniform_cortex));
                    dalpha=1.95*sqrt((total_fibers+n_fibers)/(total_fibers*n_fibers));

                    disp('Testing Uniform=Fibers Cortex')
                    if D>dalpha
                        disp(['Hypothesis rejected! With a level alpha of 0.001, percentage of diff = ' num2str((D-dalpha)*100/dalpha) '%' ])
                    else
                        disp('Hypothesis acceted!')
                    end   
                    disp(' ')

                    if show==1
                        figure(5)
                        set(figure(5),'Position',[ 1921    1    1920   1121])
                        hold on
                        plot(1:length(Emp_distr_F_cortex),Emp_distr_F_cortex,'r',1:length(Emp_distr_F_cortex),Emp_distr_Uniform_cortex,'k');
                        plot(1:length(Emp_distr_F_cortex),Emp_distr_F_cortex+dalpha,'m');
                        plot(1:length(Emp_distr_F_cortex),Emp_distr_F_cortex-dalpha,'m');
                        if D>dalpha
                            plot([I I],[0,1],'g');
                        end
                        legend('fibers','proto','accepted boundaries')
                        xlabel('X')
                        ylabel('Cumulative Probability')
                    %     bar(Density_F_B)
                    %     figurebo
                    %     bar(Density_Med_B) 
                    end

                    Write_vtk_bundles_segments(New_Points_uniform,New_number_points_curve_uniform,[],[],[],filename_fibers_uniform_down)
                else
                   disp('Density computed or parameter ''density'' set to 0')    
                end

                tElapsed = toc(tStart);
                disp(['Elapsed time: ' num2str(tElapsed) ' seconds, which means ' num2str(tElapsed/60) ' minutes' ])
                disp(' ')

                diary off

                clear('Points_Medoids','Number_points_curve_Medoids','Scalars_Medoids','Points_finale','number_points_curve_final','Scalars_color','Index_Points','Points_Cortical','Points_Basal')
            end % end BLADE
            
        end % end kind  
        fclose('all');
        lastwarn('')         
    else
        fclose('all');
        lastwarn('')
    end
end % end subject list





