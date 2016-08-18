% Reads weighted prototypes as lines in VTK format
%
% Usage: [New_Points,New_number_points_curve] = DownSampleBundlesVTK(filename,out_filename,n_fibers,data_p,data_n) 
%
% MANDATORY INPUTS:
% - filename: filename of the VTK fiber bundle to downsample
% - out_filename: filename of the downsampled  VTK fiber bundle 
% - n_fibers: number of fibers to keep 
% OPTIONAL INPUTS:
% - data_p: matrix [N,3] with the Points of the fiber bundle
% - data_n: vector [M,1] with teh number of points per streamline
%
% OUTPUTS:
% - New_Points: Matrix [P,3] with teh points of the downsampled bundle
% - New_number_points_curve: vector [n_fibers,1] with the number of points
% per streamline of the downsampled fiber bundle
%
%  Copyright Pietro GORI, Inria 
%  Written 16/08/2016

function [New_Points,New_number_points_curve] = DownSampleBundlesVTK(filename,out_filename,n_fibers,data_p,data_n) 

if nargin < 4
    disp('Downsampling - Reading Data')
    [Points,number_points_curve] = VTK_Bundles_Reader_Segments(filename);
else
    disp('Downsampling - Data already load')
    Points=data_p;
    number_points_curve=data_n;
end
            
total_fibers=length(number_points_curve);

% randomly choose n_fibers streamlines based on a uniform distribution
INDEX=randperm(total_fibers,n_fibers); 

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

Write_vtk_bundles_polyline(New_Points,New_number_points_curve,[],[],[],out_filename)

