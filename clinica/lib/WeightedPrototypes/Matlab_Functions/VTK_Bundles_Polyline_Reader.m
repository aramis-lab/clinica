% Reads weighted prototypes as lines in VTK format
%
% Usage: [Points,number_points_curve,curve_associated,Scalars] = VTK_Bundles_Polyline_Reader(filename)
%
% INPUTS:
% - filename: filename of the VTK file of the weighted prototypes saved as
% lines
%
% OUTPUTS:
% - Points: Matrix [N,3] where N is the number of points. Each row
% indicates the coordinates of a point.
% - number_points_curve: vector [M,1] where M is the number of prototypes.
% Each cell represents the number of points per prototype
% - curve_associated: vector [N,1]. Cell i indicates the prototype
% associated to the point i.
% - Scalars: vector [N,1]. Cell i indicates the weight of the prototype
% to which point i belongs to.
%
%  Copyright Pietro GORI, Inria 
%  Written 16/08/2016

function [Points,number_points_curve,curve_associated,Scalars] = VTK_Bundles_Polyline_Reader(filename)

fid = fopen(filename, 'r');
if(fid==-1)
    fprintf(1,'Error: file descriptor not valid, check the file name.\n');
    return;
end

keyWord = 'DATASET POLYDATA';
newL = GoToKeyWord(fid, keyWord);
if(newL == -1)
    fprintf(1, 'Error: file is not a vtkPolyData.\n');
    return;
end

keyWord = 'POINTS';
newL = GoToKeyWord(fid, keyWord);
if(newL==-1)
    fprintf(1, 'Cannot find flag: %s\n', keyWord);
end

% buffer = sscanf(newL,'%s%d%s');
% numPoints = buffer(length(keyWord)+1); % because these are points
numPoints = sscanf(newL,'%*s%d%*s'); % another way...
Points = zeros(numPoints,3);  % SET OF POINTS
curve_associated = zeros(numPoints,1);

newL = fgetl(fid);
count = 1;

% Read the points data
while(count<=numPoints)
    [num,c] = sscanf(newL,'%f');
    if c==3
        Points(count,:)=num';
        count=count+1;
    else
        times=c/3;
        if (mod(c,3)~=0)
            disp('Error!!')
            return
        end
        for k=0:times-1
            Points(count,:)=num(k*3+1:(k+1)*3)';
            count=count+1;            
        end
        
    end   
    newL = fgetl(fid);
end
% Points = reshape(Points, [3,numPoints])';
% end of point data

% Read the polygons
keyWord = 'LINES';
newL = GoToKeyWord(fid, keyWord);
if(newL == -1)
    return;
end

numFibers = sscanf(newL, '%*s%d%*d');
number_points_curve = zeros(numFibers,1);
F=1;

for i=1:numFibers
    newL = fgetl(fid);
    new_points = sscanf(newL,'%d');
    number_points_curve(i)=new_points(1);
    curve_associated(F:F+new_points(1)-1)=i;
    F=F+new_points(1);
end

if numPoints ~= (F-1)
    disp('Error!! numFibers is not correct')
end

% Read the scalars
Scalars = []; % set of scalars
keyWord = 'SCALARS';
newL = GoToKeyWord(fid, keyWord);
if(newL == -1)
%     fprintf(1, 'No scalar\n');
else
        
    keyWord = 'LOOKUP_TABLE';
    newL = GoToKeyWord(fid, keyWord);
    if(newL == -1)
        fprintf(1, 'No LUT\n');
        keyWord = 'SCALARS';
        newL = GoToKeyWord(fid, keyWord);
    end
    
    count = 1;
    while(count <= numPoints)
        
        newL = fgetl(fid);
        [buffer, c] = sscanf(newL, '%f');
        count = count + c;
        Scalars = [Scalars;buffer];
        
    end
    
    if (count-1)~=numPoints
        disp('Error!! count is not correct')
    end
    
end
% end of scalars

fclose(fid);