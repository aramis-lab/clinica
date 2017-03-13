function SurfStatViews( data, vol, z, layout );

%Viewer for slices of volumetric data.
%
% Usage: SurfStatViews( data, vol [, z [, layout ]] );
%
% data       = n x v matrix of data if k=1, or n x v x k array of data if 
%              k>1, or a memory map of same.
% vol.lat    = nx x ny x nz logical array, 1=in, 0=out.
% vol.vox    = 1 x 3 vector of voxel sizes in mm.
% vol.origin = position in mm of the first voxel.
% z          = 1 x s vector of z coordinates of slices. Default is 0.
%              The nearest slices are displayed. 
% layout     = matrix of indices from 1:(n*k*s) specifying which slice
%              appears in position layout(i,j), with 0 for an empty image.
%              Default is reshape(1:(n*k*s),n,k*s).

if isnumeric(data)
    [n,v,k]=size(data);
else
    datam=data;
    ss=datam.Format{2};
    if length(ss)==2
        ss=ss([2 1]);
        k=1;
    else
        ss=ss([3 1 2]);
        k=ss(3);
    end
    n=ss(1);
    v=ss(2);
end

[nx,ny,nz]=size(vol.lat);
if nargin<3 | isempty(z)
    z=0;
end
s=length(z);
if nargin<4
    layout=reshape(1:(n*k*s),n,k*s);
end
[lx,ly]=size(layout);
m=zeros(lx*ny+lx+1,ly*nx+ly+1)+NaN;
vs=[0; squeeze(cumsum(sum(sum(vol.lat))))];
for i=1:lx
    for j=1:ly
        if layout(i,j)>0
            [n1,k1,s1]=ind2sub([n,k,s],layout(i,j));
            slice=round((z(s1)-vol.origin(3))/(vol.vox(3)+(vol.vox(3)==0)))+1;
            slice=max(min(slice,nz),1);
            v1=vs(slice)+1;
            v2=vs(slice+1);
            d=zeros(nx,ny)+NaN;
            if isnumeric(data)
                d(vol.lat(:,:,slice))=data(n1,v1:v2,k1);
            else
                if length(ss)==2
                    d(vol.lat(:,:,slice))=datam.Data(1).Data(v1:v2,n1);
                else
                    d(vol.lat(:,:,slice))=datam.Data(1).Data(v1:v2,k1,n1);
                end
            end
            m((1:ny)+(ny+1)*(i-1)+1,(1:nx)+(nx+1)*(j-1)+1)=flipud(d');
        end
    end
end
m(isnan(m))=min(min(min(m(~isnan(m)))));

imagesc((1:size(m,2))*vol.vox(1),(1:size(m,1))*vol.vox(2),m); 
for i=1:lx
    for j=1:ly
        if layout(i,j)>0
            [n1,k1,s1]=ind2sub([n,k,s],layout(i,j));
            h=text(((nx+1)*(j-1)+1)*vol.vox(1),...
                ((ny+1)*(i)+1)*vol.vox(2),num2str(layout(i,j)));
            set(h,'Color','white','VerticalAlignment','bottom');
        end
    end
end
colorbar; colormap(spectral(256));
xlabel('x'); ylabel('y'); 
axis equal; axis off; 

background='white';
whitebg(gcf,background);
set(gcf,'Color',background,'InvertHardcopy','off');

set(gcf,'PaperPosition',[0.25 2.5 20 20]);%% here, you can define the image's positon and size

return
end
