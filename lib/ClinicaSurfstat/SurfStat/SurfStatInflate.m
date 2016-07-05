function surfw = SurfStatInflate( surf, w, spherefile );

%Inflates a surface towards two fitted hemi-ellipsoids.
%
% Usage: surfw = SurfStatInflate( surf [, w [, spherefile ]] );
%
% surf.coord  = 3 x v matrix of coordinates, v=#vertices.
% w           = weight in [0,1] given to hemi-ellipsoids, default 0.5.
% spherefile  = file name of a sphere surface for the left hemisphere. If
%               it is a .obj file, it assumes that the triangulation of the
%               right hemisphere is a mirror image; if it is an FS file,
%               then it assumes it is identical. The default is sphere.obj
%               for a 40962 vertex (v=81924 for both hemispheres)
%               icosahedral mesh, or lh.sphere for a 163842 vertex
%               (v=327684 for both hemispheres) icosahedral mesh.
%
% surfw.coord = 3 x v matrix of inflated coordinates.

if nargin<2
    w=0.5;
end
v=size(surf.coord,2);
if v<=81924
    if nargin<3
        spherefile='sphere.obj';
    end
    sphere=SurfStatReadSurf(spherefile);
    if v==81924
        sphere.tri=[sphere.tri sphere.tri+v];
        sphere.coord=[sphere.coord(1,:).*(sphere.coord(1,:)<0) ...
            -sphere.coord(1,:).*(sphere.coord(1,:)<0);
            sphere.coord(2:3,:) sphere.coord(2:3,:)];
    else
        if mean(surf.coord(1,:))/mean(abs(surf.coord(1,:)))<-0.5
            sphere.coord=[sphere.coord(1,:).*(sphere.coord(1,:)<0);
                sphere.coord(2:3,:)];
        else
            sphere.coord=[-sphere.coord(1,:).*(sphere.coord(1,:)<0);
                sphere.coord(2:3,:)];
        end
    end        
else
    if nargin<3
        spherefile='lh.sphere';
    end
    sphere=SurfStatReadSurf(spherefile);
    if v==327684
        sphere.tri=[sphere.tri sphere.tri+v];
        sphere.coord=[sphere.coord(1,:).*(sphere.coord(1,:)<0) ...
            sphere.coord(1,:).*(sphere.coord(1,:)>0);
            sphere.coord(2:3,:) sphere.coord(2:3,:)];
    else
        if mean(surf.coord(1,:))/mean(abs(surf.coord(1,:)))<-0.5
            sphere.coord=[sphere.coord(1,:).*(sphere.coord(1,:)<0);
                sphere.coord(2:3,:)];
        else
            sphere.coord=[sphere.coord(1,:).*(sphere.coord(1,:)>0);
                sphere.coord(2:3,:)];
        end
    end
end
maxs=max(surf.coord,[],2);
mins=min(surf.coord,[],2);
maxsp=max(sphere.coord,[],2);
minsp=min(sphere.coord,[],2);
surfw=surf;
for i=1:3
    surfw.coord(i,:)=((sphere.coord(i,:)-minsp(i))/(maxsp(i)-minsp(i))...
        *(maxs(i)-mins(i))+mins(i))*w+surf.coord(i,:)*(1-w);
end

return
end

