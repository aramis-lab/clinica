function [ surf, ab ] = SurfStatReadSurf1( filename, ab, numfields );

%Reads coordinates and triangles from a single .obj or FreeSurfer file. 
%
% Usage: [ surf, ab ] = SurfStatReadSurf1( filename [,ab [,numfields ]] );
%
% filename  = .obj or FreeSurfer file name.
% ab        = 'a' for ASCII or 'b' for binary. If it doesn't work it
%              will try the other. Default is 'a'. Ignored if FS file.
% numfields = number of fields to read, in the order below, default 4.
%
% surf.coord  = 3 x v matrix of coordinates, v=#vertices.
% surf.normal = 3 x v matrix of surface normals, only .obj files.
% surf.colr   = 4 x 1 vector of colours for the whole surface,
%               or 4 x v matrix of colours for each vertex, either 
%               uint8 in [0 255], or float in [0 1], only .obj files.
% surf.tri    = t x 3 matrix of triangle indices, 1-based, t=#triangles.
% ab          = whichever was successful.

if nargin<2
    ab='a';
end
if nargin<3
    numfields=4;
end

[pathstr,name,ext] = fileparts(filename);
if strcmp(ext,'.obj')
    % It's a .obj file
    if ab(1)=='a'
        fid=fopen(filename);
        FirstChar=fscanf(fid,'%1s',1);
        if FirstChar=='P' % ASCII
            fscanf(fid,'%f',5);
            v=fscanf(fid,'%f',1);
            surf.coord=fscanf(fid,'%f',[3,v]);
            if numfields>=2
                surf.normal=fscanf(fid,'%f',[3,v]);
                if numfields>=3
                    ntri=fscanf(fid,'%f',1);
                    ind=fscanf(fid,'%f',1);
                    if ind==0
                        surf.colr=fscanf(fid,'%f',4);
                    else
                        surf.colr=fscanf(fid,'%f',[4,v]);
                    end
                    if numfields>=4
                        fscanf(fid,'%f',ntri);
                        surf.tri=fscanf(fid,'%f',[3,ntri])'+1;
                    end
                end
            end
            fclose(fid);
        else
            fclose(fid);
            fid=fopen(filename,'r','b');
            FirstChar=fread(fid,1);
            if FirstChar==uint8(112) % binary
                fread(fid,5,'float');
                v=fread(fid,1,'int');
                surf.coord=fread(fid,[3,v],'float');
                if numfields>=2
                    surf.normal=fread(fid,[3,v],'float');
                    if numfields>=3
                        ntri=fread(fid,1,'int');
                        ind=fread(fid,1,'int');
                        if ind==0
                            surf.colr=uint8(fread(fid,4,'uint8'));
                        else
                            surf.colr=uint8(fread(fid,[4,v],'uint8'));
                        end
                        if numfields>=4
                            fread(fid,ntri,'int');
                            surf.tri=fread(fid,[3,ntri],'int')'+1;
                        end
                    end
                end
                fclose(fid);
                ab='b';
            else
                fprintf(1,'%s\n',['Unable to read ' filename ', first character ' char(FirstChar)]);
            end
        end
    else
        fid=fopen(filename,'r','b');
        FirstChar=fread(fid,1);
        if FirstChar==uint8(112) % binary
            fread(fid,5,'float');
            v=fread(fid,1,'int');
            surf.coord=fread(fid,[3,v],'float');
            if numfields>=2
                surf.normal=fread(fid,[3,v],'float');
                if numfields>=3
                    ntri=fread(fid,1,'int');
                    ind=fread(fid,1,'int');
                    if ind==0
                        surf.colr=uint8(fread(fid,4,'uint8'));
                    else
                        surf.colr=uint8(fread(fid,[4,v],'uint8'));
                    end
                    if numfields>=4
                        fread(fid,ntri,'int');
                        surf.tri=fread(fid,[3,ntri],'int')'+1;
                    end
                end
            end
            fclose(fid);
        else
            fclose(fid);
            fid=fopen(filename);
            FirstChar=fscanf(fid,'%1s',1);
            if FirstChar=='P' %ASCII
                fscanf(fid,'%f',5);
                v=fscanf(fid,'%f',1);
                surf.coord=fscanf(fid,'%f',[3,v]);
                if numfields>=2
                    surf.normal=fscanf(fid,'%f',[3,v]);
                    if numfields>=3
                        ntri=fscanf(fid,'%f',1);
                        ind=fscanf(fid,'%f',1);
                        if ind==0
                            surf.colr=fscanf(fid,'%f',4);
                        else
                            surf.colr=fscanf(fid,'%f',[4,v]);
                        end
                        if numfields>=4
                            fscanf(fid,'%f',ntri);
                            surf.tri=fscanf(fid,'%f',[3,ntri])'+1;
                        end
                    end
                end
                fclose(fid);
                ab='a';
            else
                fprintf(1,'%s\n',['Unable to read ' filename ', first character ' char(FirstChar)]);
            end
        end
    end
else
    % Assume it's a FreeSurfer file
    fid = fopen(filename, 'rb', 'b') ;
    b1 = fread(fid, 1, 'uchar') ;
    b2 = fread(fid, 1, 'uchar') ;
    b3 = fread(fid, 1, 'uchar') ;
    magic = bitshift(b1, 16) + bitshift(b2,8) + b3 ;
    if magic==16777214
        fgets(fid);
        fgets(fid);% read the first two lines
        v = fread(fid, 1, 'int32') ; % number of vertex %FS pial, 163842
        t = fread(fid, 1, 'int32') ;% number of triangles %FS pial, 327680
        surf.coord = fread(fid, [3 v], 'float32') ; % the coordinate for all the vertex
        if numfields==4
            surf.tri = fread(fid, [3 t], 'int32')' + 1 ;
        end
        fclose(fid) ;
    else
        fprintf(1,'%s\n',['Unable to read ' filename ', magic = ' num2str(magic)]);
    end
    ab='b';
end

return
end

