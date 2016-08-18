

function Write_vtk_bundles_polyline(Points,number_points_curve,Scalars, Colors, TextureCoordinates,filename)

[R,C] = size(Points);

if C ~= 3
    disp('Error, Points must be a matrix [R,3]')
    return;
end

number_points = R;
number_fibers = length(number_points_curve);

fid = fopen(filename, 'w');

fprintf(fid, '# vtk DataFile Version 3.0\nvtk output\nASCII\n');
fprintf(fid, 'DATASET POLYDATA\n');
fprintf(fid, 'POINTS %d float\n',number_points);

for i=1:number_points
    fprintf(fid,'%10.20f %10.20f %10.20f\n',Points(i,1),Points(i,2),Points(i,3));
end

fprintf(fid,'\nLINES %d %d\n ', number_fibers, number_points+number_fibers);
k=0;

for i=1:length(number_points_curve)   
    fprintf(fid,'%d ',number_points_curve(i));
    for j=1:number_points_curve(i)-1        
        fprintf(fid,'%d ',k);
        k=k+1;                   
    end
    fprintf(fid,'%d\n ',k); 
    k=k+1;
end

fprintf(fid,'\nPOINT_DATA %d \n',number_points);

if( ~isempty(Scalars) ) % Scalars is a vector [nPts,1], to each point it gives a color

        if (length(Scalars)~=number_points)
            disp('Warning: not point data!');
        end

    fprintf(fid,'SCALARS scalars double\n'); % it specifies that we use double and that it is called scalars
    fprintf(fid, 'LOOKUP_TABLE default\n'); % we use the default LOOKUP_TABLE and each scalar says the colour assigned to each point
    for i=1:length(Scalars)
        fprintf(fid, '%10.20f\n',Scalars(i));
    end      
end

if (~isempty(Colors)) % write colors if any

    if (size(Colors,1)~=number_points)
        disp('Warning: nb of columns in colors array does not match number of points!');
    end

    nValues = size(Colors,2);
    
    fprintf(fid,'COLOR_SCALARS colors %d\n', nValues);
    for i=1:number_points
        for j=1:nValues
            fprintf(fid,'%d ', Colors(i,j));
        end
        fprintf(fid,'\n');
    end

end

if (~isempty(TextureCoordinates)) 

    if (size(TextureCoordinates,1)~=number_points)
        disp('Warning: nb of columns in TextureCoordinates array does not match number of points!');
    end

    nValues = size(TextureCoordinates,2);
   
    fprintf(fid,'TEXTURE_COORDINATES TextureCoordinates %d float\n', nValues);
    for i=1:number_points
        for j=1:nValues
            fprintf(fid,'%f ', TextureCoordinates(i,j));
        end
        fprintf(fid,'\n');
    end

end

fclose(fid);

end