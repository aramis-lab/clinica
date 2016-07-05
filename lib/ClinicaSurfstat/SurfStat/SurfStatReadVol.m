function [data,vol] = SurfStatReadVol(filenames,mask,step,dirname,maxmem);

%Reads volumetric files in MINC, ANALYZE, NIFTI or AFNI format. 
% 
% Usage: [ data, vol ] = SurfStatReadVol( filenames [,mask [,step  ...
%                        [,dirname [, maxmem ] ] ] ] ). 
%
% filenames = file name with extension .mnc, .img, .nii or .brik as above 
%             (n=1), or n x k cell array of file names. Either the data
%             files are 4D, or k>1, but not both. If the data files are 4D
%             then k=size of the 4th dimension.  
% mask      = 3D logical array of the same size as the unsampled data,
%             1=inside, 0=outside. Default is ones, i.e. all voxels.
% step      = 1 x 3 vector of steps for subsampling the volume; only data
%             at indices 1, 1+step, 1+2*step, .... is stored. If a scalar,
%           = [step step step], or 3 x 3 matrix of 
%           = [startx stepx stopx; starty stepy stopy; startz stepz stopz]. 
% dirname   = name of a directory where you have write permission. This
%             is only needed if the data is memory mapped. The data is
%             memory mapped only if it exceeds maxmem Mb as 4 byte reals.
%             SurfStatReadVol then writes the data as 4 byte reals to a
%             temporary file in dirname, then memory maps this file to the
%             output data. If dirname does not exist, SurfStatReadVol will
%             create it. The default is
%           = (filenames{1,1} directory)/SurfStat, so you can ignore this
%             parameter if you have write permission in the filenames{1,1} 
%             directory.
% maxmem    = memory limit in Mb. The default is 64. 
%
% If filenames is a cell array of file names then 
% data       = n x v matrix of data if k=1, or n x v x k array of data if 
%              k>1, or a memory map of same. Note that the mapped file is
%              not deleted after you quit MATLAB.
% vol.lat    = 3D logical array, 1=in, 0=out, clamped to the mask.
% vol.vox    = 1 x 3 vector of voxel sizes in mm of the clamped mask.
% vol.origin = position in mm of the first voxel of the clamped mask.
% vol.parent = header info from the first file.
% If n=1, data is double precision; 
% if n>1, data is single precision.

if isstr(filenames)
    filenames={filenames};
end

global MAXMEM
if nargin<5
    if isempty(MAXMEM)
        maxmem=64;
    else
        maxmem=MAXMEM;
    end
end
maxmem=maxmem*2^20;

d=SurfStatReadVol1(filenames{1,1},0,0);
n=size(filenames,1);
if length(d.dim)==3 | d.dim(4)==1
    d.dim(4)=1;
    k=size(filenames,2);
else
    k=d.dim(4);
    kk=size(filenames,2);
    if kk>1 & k>1
        error(['Can''t read k=' num2str(kk) '>1 files of 4D data, dims= ' num2str(d.dim)]);
    end
end

if nargin<2 | isempty(mask)
    mask=logical(ones(d.dim(1:3)));
end
if nargin<3 | isempty(step)
    step=1;
end
if length(step)==1
    step=[step step step];
end
if isnumeric(step)
    i0=1:step(1):d.dim(1);
    j0=1:step(2):d.dim(2);
    k0=1:step(3):d.dim(3);
else
    i0=round((step{1}-d.origin(1))/d.vox(1)+1);
    j0=round((step{2}-d.origin(2))/d.vox(2)+1);
    k0=round((step{3}-d.origin(3))/d.vox(3)+1);
    if isempty(i0)
        i0=1:d.dim(1);
    else
        i0=max(min(i0,d.dim(1)),1);
    end
    if isempty(j0)
        j0=1:d.dim(2);
    else
        j0=max(min(j0,d.dim(2)),1);
    end
    if isempty(k0)
        k0=1:d.dim(3);
    else
        k0=max(min(k0,d.dim(3)),1);
    end
end    
vol.lat=mask(i0,j0,k0);
mi=find(any(any(vol.lat,2),3));
mj=find(any(any(vol.lat,1),3));
mk=find(any(any(vol.lat,1),2));
bbox=[min(mi) max(mi); min(mj) max(mj); min(mk) max(mk)]';
vol.lat=vol.lat(bbox(1,1):bbox(2,1),bbox(1,2):bbox(2,2),bbox(1,3):bbox(2,3));
steps=zeros(1,3);
if length(i0)>1
    steps(1)=i0(2)-i0(1);
end
if length(j0)>1
    steps(2)=j0(2)-j0(1);
end
if length(k0)>1
    steps(3)=k0(2)-k0(1);
end
vol.vox=d.vox(1:3).*steps;
vol.origin=d.origin+(bbox(1,:)-1).*vol.vox;
vol.parent=d;

slices=k0(1)+(mk-1)*steps(3);
slicedlat=logical(zeros(d.dim(1),d.dim(2),length(slices)));
is=i0(1)+((bbox(1,1):bbox(2,1))-1).*steps(1);
js=j0(1)+((bbox(1,2):bbox(2,2))-1).*steps(2);
ks=mk-bbox(1,3)+1;
slicedlat(is,js,:)=vol.lat(:,:,ks);

if n>1
    fprintf(1,'%s',[num2str(n) ' files to read, % remaining: 100 ']);
    n10=floor(n/10);
end
v=sum(vol.lat(:));
isnum=(n*v*k*4<=maxmem);
if isnum
    data=zeros(n,v,k,'single');
else
    if nargin<4 | isempty(dirname)
        [PATHSTR,NAME,EXT]=fileparts(filenames{1});
        dirname=fullfile(PATHSTR,'SurfStat');
    end
    if ~exist(dirname,'dir')
        [SUCCESS,MESSAGE,MESSAGEID]=mkdir(dirname);
        if ~SUCCESS
            error(MESSAGEID,['Tried to make directory ' dirname ' for memory mapping. \n',...
                'Please specify the name of a directory where you have write permission \n',...
                'as the fourth parameter of SurfStatReadVol - see the help :-) Error: \n',...
                MESSAGE]);
        end
    end
    [PATHSTR,NAME,EXT]=fileparts(tempname);
    Filename=fullfile(dirname,NAME);
    fid=fopen(Filename,'wb');
end
for i=1:n
    if n>1 & rem(i,n10)==0
        fprintf(1,'%s',[num2str(round(100*(1-i/n))) ' ']);
    end
    for j=1:k
        if d.dim(4)==1
            d=SurfStatReadVol1(filenames{i,j},slices,1);
        else
            d=SurfStatReadVol1(filenames{i},slices,j);
        end
        data2=d.data(slicedlat)';
        data2(isnan(data2))=0;
        if isnum
            data(i,:,j)=data2;
        else
            fwrite(fid,data2,'single');
        end
    end
end
if n>1
    fprintf(1,'%s\n','Done');
end
if ~isnum
    fclose(fid);
    if k==1
        data=memmapfile(Filename,'Format',{'single' [v n] 'Data'},'Writable',true)
    else
        data=memmapfile(Filename,'Format',{'single' [v k n] 'Data'},'Writable',true)
    end
end

return;
end

