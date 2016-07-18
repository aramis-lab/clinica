function SurfStatWriteSurf( filenames, surf, ab );

%Writes coordinates and triangles to an array of .obj or FreeSurfer files. 
%
% Usage: SurfStatWriteSurf( filenames, surf [,ab] );
% 
% filenames = .obj or FS file name (n=1) or n x k cell array of file names.
%             If extension=.obj, writes a .obj file (ASCII or binary), 
%             else writes a FS file (binary only).
%             If k>1 then the data is split equally between the k files.
% surf.coord = n x v x 3 array of coordinates, v=#vertices, 3 x v if n=1.
% surf.tri   = t x 3 matrix of triangle indices, 1-based, t=#triangles.
% ab         = 'a' for ASCII (default) or 'b' for binary.
%
% For .obj files, the following are optional:
% surf.normal = n x v x 3 array of surface normals, 3 x v if n=1.
% surf.colr   = n x 4 matrix of colours for the whole surface, 
%            or n x v x 4 aray of colours for each vertex, 4 x v if n=1,
%        either uint8  in [0 255], or float in [0 1]. Default is ones(n,4).

if nargin<3
    ab='a';
end

if isstr(filenames)
    sf=filenames;
    filenames=cell(1,1);
    filenames(1)={sf};
end
[n,k]=size(filenames);

v=size(surf.coord,2+(n>1));
t=size(surf.tri,1);

if rem(v,k)>0 | rem(t,k)>0 
    return
end

vk=round(v/k);
v1=1:vk:v;
v2=vk:vk:v;
tk=round(t/k);
t1=1:tk:t;
t2=tk:tk:t;

for j=1:k
    surf1.tri=surf.tri(t1(j):t2(j),:)-v1(j)+1;
    for i=1:n
        if n==1
            surf1.coord=surf.coord(:,v1(j):v2(j));
            if isfield(surf,'colr')
                if size(colr,1)==1
                    surf1.colr=surf.colr';
                else
                    surf1.colr=surf.colr(:,v1(j):v2(j));
                end
            end
            if isfield(surf,'normal')
                surf1.normal=surf.normal(:,v1(j):v2(j));
            end
        else
            surf1.coord=squeeze(surf.coord(i,v1(j):v2(j),:))';
            if isfield(surf,'colr')
                if ndims(colr)==2
                    surf1.colr=surf.colr(i,:)';
                else
                    surf1.colr=squeeze(surf.colr(i,v1(j):v2(j),:))';
                end
            end
            if isfield(surf,'normal')
                surf1.normal=squeeze(surf.normal(i,v1(j):v2(j),:))';
            end
        end
        SurfStatWriteSurf1(filenames{i,j},surf1,ab);
    end
end

return
end


