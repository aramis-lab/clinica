function [New_Points,New_number_points_curve] = DownSampleBundlesVTK_uniform_number_fast(filename,new_filename,n_fibers,data_p,data_n) 

if nargin < 4
    disp('Downsampling - Reading Data')
    [Points,number_points_curve] = VTK_Bundles_Reader_Segments(filename);
else
    disp('Downsampling - Data already load')
    Points=data_p;
    number_points_curve=data_n;
end
            
total_fibers=length(number_points_curve);

INDEX=randperm(total_fibers,n_fibers); % faster...

New_number_points_curve=number_points_curve(INDEX);
Total_New_Points=sum(New_number_points_curve);

New_Points = zeros(Total_New_Points,3);
l=1;

for j=1:length(INDEX)
    k=1+sum(number_points_curve(1:(INDEX(j)-1)));
    number_points=number_points_curve(INDEX(j));
    New_Points(l:l+number_points-1,:)=Points(k:k+number_points-1,:);
    l=l+number_points;
end

%Write_vtk_bundles_segments(New_Points,New_number_points_curve,[],[],[],[new_filename(1:end-4) '_seg_' num2str(n_fibers) '.vtk'])
Write_vtk_bundles_polyline(New_Points,New_number_points_curve,[],[],[],[new_filename(1:end-4) '_poly_' num2str(n_fibers) '.vtk'])

