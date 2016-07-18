function s = SurfStatVol2Surf( vol, surf );

%Interpolates a volume to each surface, then averages. 
%
% Usage: s = SurfStatVol2Surf( vol, surfs );
%
% vol.data   = nx x ny x nz volume aligned to the surface data.
% vol.origin = [x,y,z] location of vol.data(1,1,1) in mm.
% vol.vox    = [x,y,z] voxel size in mm.
% surf.coord = 3 x v matrix of coordinates in mm, v=#vertices, in double
%              precision, if n=1; or n x v x 3 in single precision if n>1,
%              or memory map of same.
%
% s = 1 x v vector of vol.data interpolated to each surface, then averaged.

if isnumeric(surf.coord)
    v=size(surf.coord,2);
    one=ones(v,1);
    if ndims(surf.coord)==2
        vox=(double(surf.coord')-one*vol.origin)./(one*vol.vox(1:3))+1;
        s=interpn(vol.data,vox(:,1),vox(:,2),vox(:,3),'linear',0);
    else
        n=size(surf.coord,1);
        s=zeros(v,1);
        fprintf(1,'%s',[num2str(n) ' surfaces to interpolate, % remaining: 100 ']);
        n10=floor(n/10);
        for i=1:n
            if rem(i,n10)==0
                fprintf(1,'%s',[num2str(round(100*(1-i/n))) ' ']);
            end
            c=double(squeeze(surf.coord(i,:,:)));
            vox=(c-one*vol.origin)./(one*vol.vox(1:3))+1;
            s=s+interpn(vol.data,vox(:,1),vox(:,2),vox(:,3),'linear',0);
        end
        s=s'/n;
        fprintf(1,'%s\n','Done');
    end
else
    sz=surf.coord.Format{2};
    v=sz(1);
    n=sz(3);
    one=ones(v,1);
    s=zeros(v,1);
    fprintf(1,'%s',[num2str(n) ' surfaces to interpolate, % remaining: 100 ']);
    n10=floor(n/10);
    for i=1:n
        if rem(i,n10)==0
            fprintf(1,'%s',[num2str(round(100*(1-i/n))) ' ']);
        end
        c=double(surf.coord.Data(1).Data(:,:,i));
        vox=(c-one*vol.origin)./(one*vol.vox(1:3))+1;
        s=s+interpn(vol.data,vox(:,1),vox(:,2),vox(:,3),'linear',0);
    end
    s=s'/n;
    fprintf(1,'%s\n','Done');
end

return
end
