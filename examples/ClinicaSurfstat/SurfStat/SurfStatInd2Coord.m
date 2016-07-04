function coord = SurfStatInd2Coord( ind, surf );

%Converts a vertex index to x,y,z coordinates.
%
% Usage: coord = SurfStatInd2Coord( ind, surf );
%
% ind          = 1 x c vector of indices of vertex, 1-based.
% surf.coord   = 3 x v matrix of coordinates.
% or
% surf.lat     = 3D logical array, 1=in, 0=out. 
% surf.vox     = 1 x 3 vector of voxel sizes in mm, [1 1 1] by default.
% surf.origin  = position in mm of the first voxel, [0 0 0] by default.
%
% coord = 3 x c matrix of coordinates. 

if isfield(surf,'coord')
    coord=surf.coord(:,ind);
end
if isfield(surf,'lat')
    if ~isfield(surf,'vox')
        surf.vox=ones(1,3);
    end
    if ~isfield(surf,'origin');
        surf.origin=zeros(1,3);
    end
    vid=cumsum(surf.lat(:)).*surf.lat(:);
    [tf,loc]=ismember(ind,vid);
    loc(loc==0)=NaN;
    dim=size(surf.lat);
    [i,j,k]=ind2sub(dim,loc);
    coord=zeros(3,length(ind));
    coord(1,:)=surf.origin(1)+(i-1)*surf.vox(1);
    coord(2,:)=surf.origin(2)+(j-1)*surf.vox(2);
    coord(3,:)=surf.origin(3)+(k-1)*surf.vox(3);
end

return
end

