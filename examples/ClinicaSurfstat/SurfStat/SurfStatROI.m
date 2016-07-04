function maskROI = SurfStatROI( centre, radius, surf );

%ROI on the surface or volume with a given centre and radius.
%
% Usage: ROImask = SurfStatROI( centre, radius, surf );
% 
% centre      = id number of vertex, or 
%               3 x 1 vector of [x; y; z] coordinates in mm.
% radius      = radius, mm.
% surf.coord  = 3 x v matrix of coordinates of surface.
% or
% surf.lat    = nx x ny x nz array, 1=in, 0=out, clamped to the mask.
% surf.vox    = 1 x 3 vector of voxel sizes in mm of the clamped mask.
% surf.origin = position in mm of the first voxel of the clamped mask.
%
% maskROI = 1 x v vector, 1=inside ROI, 0=outside.

if isfield(surf,'coord')
    if length(centre)==1
        id=centre;
        centre=surf.coord(:,id);
    end
    d2=sum((centre*ones(1,size(surf.coord,2))-surf.coord).^2);
    maskROI=d2<radius^2;
else
    if length(centre)==1
        id=centre;
        vid=int32(cumsum(surf.lat(:)).*surf.lat(:));
        [i,j,k]=ind2sub(find(vid==id),size(surf.lat));
        centre=(([i,j,k]-1).*surf.vox+surf.origin)';
    end
    dim=size(surf.lat);
    x0=((1:dim(1))-1)*surf.vox(1)+surf.origin(1);
    y0=((1:dim(2))-1)*surf.vox(2)+surf.origin(2);
    z0=((1:dim(3))-1)*surf.vox(3)+surf.origin(3);
    [i,j,k]=ndgrid(x0,y0,z0);
    coord=[i(surf.lat) j(surf.lat) k(surf.lat)]';
    d2=sum((centre*ones(1,size(coord,2))-coord).^2);
    maskROI=d2<radius^2;
end

return
end
    