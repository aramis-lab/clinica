function ind = SurfStatCoord2Ind( coord, surf );

%Converts a vertex index to x,y,z coordinates.
%
% Usage: coord = SurfStatInd2Coord( ind, surf );
%
% coord        = c x 3 matrix of coordinates for finding indices thereof. 
% surf.coord   = 3 x v matrix of coordinates.
% or
% surf.lat     = 3D logical array, 1=in, 0=out. 
% surf.vox     = 1 x 3 vector of voxel sizes in mm, [1 1 1] by default.
% surf.origin  = position in mm of the first voxel, [0 0 0] by default.
%
% ind = c x 1 vector of indices of the nearest vertex to the surface, 
%       1-based. If surf is a volume and the point is outside, then ind=0.

c=size(coord,1);
ind=zeros(c,1);
if isfield(surf,'coord')
    v=size(surf.coord,2);
    for i=1:c
       dist=sum((surf.coord-repmat(coord(i,:)',1,v)).^2);
       ind(i)=find(dist==min(dist));
    end
end
if isfield(surf,'lat')
    if ~isfield(surf,'vox')
        surf.vox=ones(1,3);
    end
    if ~isfield(surf,'origin');
        surf.origin=zeros(1,3);
    end
    i=round((coord(:,1)-surf.origin(1))/(surf.vox(1)+(surf.vox(1)==0))+1);
    j=round((coord(:,2)-surf.origin(2))/(surf.vox(2)+(surf.vox(2)==0))+1);
    k=round((coord(:,3)-surf.origin(3))/(surf.vox(3)+(surf.vox(3)==0))+1);
    dim=size(surf.lat);
    i(i<1 | i>dim(1))=0;
    j(j<1 | j>dim(2))=0;
    k(k<1 | k>dim(3))=0;
    a=i&j&k;
    ind=zeros(c,1);
    vid=cumsum(surf.lat(:)).*surf.lat(:);
    ind(a)=vid(sub2ind(dim,i(a),j(a),k(a)));
end

return
end

