function vol = SurfStatSurf2Vol( s, surf, template );

%Interpolates data on each surface to a volume, then averages. 
%
% Usage: vol = SurfStatSurf2Vol( s, surf [, template] );
%
% s          = 1 x v vector of surface data, v=#vertices.
% surf.coord = 3 x v matrix of coordinates in mm, v=#vertices, in double
%              precision, if n=1; or n x v x 3 in single precision if n>1,
%              or memory map of same.
% template   = template volume file, icbm_template_2.00mm.mnc by default. 
%              For ANALYZE format use icbm_template_2.00mm.img,
%              for NIFTI   format use icbm_template_2.00mm.nii,
%              for AFNI    format use icbm_template_2.00mm.brik.
%
% vol.data   = nx x ny x nz array; for each surf.coord, s is assigned to 
%              the nearest voxel in vol.data, then averaged over surfaces.
% vol.origin = [x,y,z] location of vol.data(1,1,1) in mm.
% vol.vox    = [x,y,z] voxel size in mm.

if nargin<3
    template='icbm_template_2.00mm.mnc';
end

vol=SurfStatReadVol1(template);
[nx,ny,nz]=size(vol.data);

if isnumeric(surf.coord)
    v=size(surf.coord,2);
    one=ones(v,1);
    if ndims(surf.coord)==2
        vox=round((double(surf.coord')-one*vol.origin)./(one*vol.vox(1:3))+1);
        in=vox(:,1)>=1&vox(:,2)>=1&vox(:,3)>=1;
        in=in&vox(:,1)<=nx&vox(:,2)<=ny&vox(:,3)<=nz;
        vol.data=accumarray(vox(in,:),s(in),[nx,ny,nz]);
        ns=accumarray(vox(in,:),1,[nx,ny,nz]);
        vol.data=vol.data./(ns+(ns<=0)).*(ns>0);
    else
        n=size(surf.coord,1);
        ns=zeros([nx,ny,nz]);
        fprintf(1,'%s',[num2str(n) ' surfaces to interpolate, % remaining: 100 ']);
        n10=floor(n/10);
        for i=1:n
            if rem(i,n10)==0
                fprintf(1,'%s',[num2str(round(100*(1-i/n))) ' ']);
            end
            c=double(squeeze(surf.coord(i,:,:)));
            vox=round((c-one*vol.origin)./(one*vol.vox(1:3))+1);
            in=vox(:,1)>=1&vox(:,2)>=1&vox(:,3)>=1;
            in=in&vox(:,1)<=nx&vox(:,2)<=ny&vox(:,3)<=nz;
            vs=accumarray(vox(in,:),s(in),[nx,ny,nz]);
            ns=accumarray(vox(in,:),1,[nx,ny,nz]);
            vol.data=vol.data+vs./(ns+(ns<=0)).*(ns>0);
        end
        vol.data=vol.data/n;
        fprintf(1,'%s\n','Done');
    end
    vol.file_name=inputname(1);
else
    sz=surf.coord.Format{2};
    v=sz(1);
    n=sz(3);
    one=ones(v,1);
    ns=zeros([nx,ny,nz]);
    fprintf(1,'%s',[num2str(n) ' surfaces to interpolate, % remaining: 100 ']);
    n10=floor(n/10);
    for i=1:n
        if rem(i,n10)==0
            fprintf(1,'%s',[num2str(round(100*(1-i/n))) ' ']);
        end
        c=double(surf.coord.Data(1).Data(:,:,i));
        vox=round((c-one*vol.origin)./(one*vol.vox(1:3))+1);
        in=vox(:,1)>=1&vox(:,2)>=1&vox(:,3)>=1;
        in=in&vox(:,1)<=nx&vox(:,2)<=ny&vox(:,3)<=nz;
        vs=accumarray(vox(in,:),s(in),[nx,ny,nz]);
        ns=accumarray(vox(in,:),1,[nx,ny,nz]);
        vol.data=vol.data+vs./(ns+(ns<=0)).*(ns>0);
    end
    vol.data=vol.data/n;
    fprintf(1,'%s\n','Done');
    vol.file_name=inputname(1);
end    
    
return
end
