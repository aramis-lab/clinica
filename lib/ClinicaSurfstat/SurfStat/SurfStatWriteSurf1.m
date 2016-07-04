function surf = SurfStatWriteSurf1( filename, surf, ab );

%Writes coordinates and triangles to a single .obj or FreeSurfer file.
%
% Usage: surf = SurfStatWriteSurf1( filename, surf [,ab] );
%
% filename   = .obj or FS file name. If extension=.obj, writes a .obj file
%              (ASCII or binary), else writes a FS file (binary only).
% surf.coord = 3 x v matrix of coordinates, v=#vertices.
% surf.tri   = t x 3 matrix of triangle indices, 1-based, t=#triangles.
% ab         = 'a' for ASCII (default) or 'b' for binary.
%
% For .obj files, the following are optional and returned if not present:
% surf.normal = 3 x v matrix of surface normals.
% surf.colr   = 4 x 1 vector of colours for the whole surface,
%            or 4 x v matrix of colours for each vertex, either uint8
%            in [0 255], or float in [0 1]. Default is [1 1 1 1]'.

v=size(surf.coord,2);
t=size(surf.tri,1);

[pathstr,name,ext] = fileparts(filename);
if strcmp(ext,'.obj')
    % It's a .obj file
    
    if ~isfield(surf,'normal')
        u1=surf.coord(:,surf.tri(:,1));
        d1=surf.coord(:,surf.tri(:,2))-u1;
        d2=surf.coord(:,surf.tri(:,3))-u1;
        c=cross(d1,d2,1);
        surf.normal=zeros(3,v);
        for j=1:3
            for k=1:3
                surf.normal(k,:)=surf.normal(k,:)+accumarray(surf.tri(:,j),c(k,:)')';
            end
        end
        surf.normal=surf.normal./(ones(3,1)*sqrt(sum(surf.normal.^2,1)));
    end
    if ~isfield(surf,'colr')
        surf.colr=[1 1 1 1]';
    end

    if nargin<3 | ab(1)=='a'
        fid=fopen(filename,'w');
        fprintf(fid,'P 0.3 0.3 0.4 10 1 %d \n',v);
        fprintf(fid,'%f %f %f \n',surf.coord);
        fprintf(fid,'  \n');
        fprintf(fid,'%f %f %f \n',surf.normal);
        fprintf(fid,'  \n');
        fprintf(fid,'%d %d \n',t,(size(surf.colr,2)>1)*2);
        if isa(surf.colr,'uint8')
            fprintf(fid,'%f %f %f %f \n',double(surf.colr)/255);
        else
            fprintf(fid,'%f %f %f %f \n',surf.colr);
        end
        fprintf(fid,'  \n');
        fprintf(fid,'%d %d %d %d %d %d %d %d \n',(1:t)*3);
        fprintf(fid,'  \n');
        fprintf(fid,'%d %d %d %d %d %d %d %d \n',surf.tri'-1);
        fclose(fid);
    else
        fid=fopen(filename,'w','b');
        fwrite(fid,uint8(112),'uint8');
        fwrite(fid,[0.3 0.3 0.4 10 1],'float');
        fwrite(fid,v,'int');
        fwrite(fid,surf.coord,'float');
        fwrite(fid,surf.normal,'float');
        fwrite(fid,t,'int');
        fwrite(fid,(size(surf.colr,2)>1)*2,'int');
        if isa(surf.colr,'uint8')
            fwrite(fid,surf.colr,'uint8');
        else
            fwrite(fid,uint8(round(surf.colr*255)),'uint8');
        end
        fwrite(fid,(1:t)*3,'int');
        fwrite(fid,surf.tri'-1,'int');
        fclose(fid);
    end
else
    % Assume it's a FreeSurfer file
    fid = fopen(filename, 'wb', 'b') ;
    magic = 16777214;
    b1 = bitand(bitshift(magic, -16), 255) ;
    b2 = bitand(bitshift(magic, -8), 255) ;
    b3 = bitand(magic, 255) ;
    fwrite(fid, [b1 b2 b3], 'uchar') ;
    fwrite(fid, ['Created by SurfStat on ' datestr(now) char(10) char(10)], 'char');
    fwrite(fid, [v t], 'int32') ;
    fwrite(fid, surf.coord, 'float') ;
    fwrite(fid, surf.tri'-1, 'int32') ;
    fclose(fid);
end

return


