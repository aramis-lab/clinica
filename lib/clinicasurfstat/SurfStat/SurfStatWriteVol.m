function SurfStatWriteVol( filenames, data, vol );

%wtites volumetric data to files in MINC, ANALYZE, NIFTI or AFNI format. 
%
% Usage: SurfStatWriteVol( filenames, data, vol );
%
% filenames  = single file name with extension .mnc, .img, .nii or .brik as 
%              above (3D if k=1, 4D if k>1), or cell array of k 3D files. 
% data       = k x v matrix of data, v=#vertices.
% vol.lat    = nx x ny x nz array, 1=in, 0=out, clamped to the mask.
% vol.vox    = 1 x 3 vector of voxel sizes in mm of the clamped mask.
% vol.origin = position in mm of the first voxel of the clamped mask.
% vol.parent = header info from SurfStatReadVol1.

n=size(vol.lat);
m=zeros(n);
if ~isfield(vol,'parent_file')
    d=vol.parent;
end
d.data=zeros(d.dim(1:3));
step=round(vol.vox./d.vox);
s=round((vol.origin-d.origin)./d.vox)+1;
[k,v]=size(data);
for i=1:k
    m(vol.lat)=data(i,:);
    d.data(s(1)+(0:n(1)-1)*step(1),s(2)+(0:n(2)-1)*step(2),...
        s(3)+(0:n(3)-1)*step(3))=m;
    if isstr(filenames)
        d.file_name=filenames;
        SurfStatWriteVol1(d,1:d.dim(3),i);
    else
        d.file_name=filenames{i};
        SurfStatWriteVol1(d,1:d.dim(3),1);
    end
end

return
end






