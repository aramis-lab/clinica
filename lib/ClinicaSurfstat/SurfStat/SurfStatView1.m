function SurfStatView1(struct, surf, varargin );

%Viewer for surface and/or volume data or P-value or Q-value structures.
% 
% Usage: SurfStatView( struct, surf, ... );
% 
% struct        = 1 x v vector of data, v=#vertices,
%                 or one of the following structures:
% P-value structure:
% struct.P      = 1 x v vector of corrected P-values for vertices.
% struct.C      = 1 x v vector of corrected P-values for clusters.
% struct.mask   = 1 x v, 1=inside, 0=outside, v=#vertices.
% struct.thresh = P-value threshold for plot, 0.05 by default.
% Q-value structure:
% struct.Q      = 1 x v matrix of Q-values. 
% struct.mask   = 1 x v, 1=inside, 0=outside, v=#vertices.
% struct.thresh = Q-value threshold for plot, 0.05 by default.
%
% surf is either a surface or volume structure:
% Surface structure:
% surf.coord   = 3 x v matrix of surface coordinates.
% surf.tri     = t x 3 matrix of triangle indices, 1-based, t=#triangles.
% Volume structure:
% surf.lat     = 3D logical array, 1=in, 0=out.
% surf.vox     = 1 x 3 vector of voxel sizes in mm, [1 1 1] by default.
% surf.origin  = position in mm of the first voxel, [0 0 0] by default.
% 
% The following can be added as extra pairs of parameters in any order: 
%
% For slices of volume structures:
% 'x', vector of slice positions in x direction, default mean(x coord).
% 'y', vector of slice positions in y direction, default mean(y coord).
% 'z', vector of slice positions in z direction, default mean(z coord).
%
% For isosurfaces of volume structures: 
% 'maskthresh' mask threshold, default 0.5.
% 'datathresh' data threshold, default [].
% 'maskalpha', mask    alpha, default 0.2.
% 'dataalpha', data    alpha, default 1.
% 'clusalpha', cluster alpha, default 0.2.
% 'maskcolor', mask    colour, any matlab ColorSpec, default 'green'.
% 'datacolor', data    colour, any matlab ColorSpec, default 'red'.
% 'cluscolor', cluster colour, any matlab ColorSpec, default 'blue'.
%
% 'title', any string, data name by default.
% 'background', background colour, any matlab ColorSpec, such as 
%   'white' (default), 'black'=='k', 'r'==[1 0 0], [1 0.4 0.6] (pink) etc.
%   Letter and line colours are inverted if background is dark (mean<0.5). 
%
% E.g. SurfStatView( struct, surf, 'x', 10, 'datathresh', 3, ...
%          'title', 'One x slice at 10mm, thresholded at 3' );
%
% To change the colour limits to e.g. [2 4], use set( gca, 'CLim', [2 4] );
% To change the colour map, use e.g. colormap( 'jet' ). 
% Surfaces can be edited in the figure window by clicking e.g. "Rotate 3D".

mask=[];
maskthresh=0.5;
datathresh=[];
maskalpha=0.2;
dataalpha=1;
clusalpha=0.2;
maskcolor='green';
datacolor='red';
cluscolor='blue';
xs=NaN;
ys=NaN;
zs=NaN;
title1=inputname(1);
background='white';

optargin=size(varargin,2);
for narg=1:2:optargin
    switch varargin{narg}
        case 'mask'
            mask=varargin{narg+1};
        case 'maskthresh'
            maskthresh=varargin{narg+1};
        case 'datathresh'
            datathresh=varargin{narg+1};
        case 'maskalpha'
            maskalpha=varargin{narg+1};
        case 'dataalpha'
            dataalpha=varargin{narg+1};
        case 'clusalpha'
            clusalpha=varargin{narg+1};
        case 'maskcolor'
            maskcolor=varargin{narg+1};
        case 'datacolor'
            datacolor=varargin{narg+1};
        case 'cluscolor'
            cluscolor=varargin{narg+1};
        case 'x'
            xs=varargin{narg+1};
        case 'y'
            ys=varargin{narg+1};
        case 'z'
            zs=varargin{narg+1};
        case 'title'
            title1=varargin{narg+1};
        case 'background'
            background=varargin{narg+1};
        otherwise
            warning(['Unrecognized optional argument ' varargin{narg} ' ignored :-(']);
    end
end

hold on;

if isfield(surf,'coord')
    if isempty(struct)
        struct=zeros(1,size(surf.coord,2));
    end
    
    if ~isstruct(struct)
        clim=[min(struct),max(struct)];
        if clim(1)==clim(2)
            clim=clim(1)+[-1 0];
        end
        a=trisurf(surf.tri,surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),...
            double(struct),'EdgeColor','none','FaceColor','interp');  
        colorbar; colormap(spectral(256));
        set(gca,'CLim',clim);
        set(gca,'Tag','SurfStatView 1 0');
        dcm_obj=datacursormode(gcf);
        set(dcm_obj,'UpdateFcn',@SurfStatDataCursor,'DisplayStyle','window');
    else
        if isfield(struct,'P')
            if ~isfield(struct,'thresh')
                struct.thresh=0.05;
            end

            signifpeak=struct.P<struct.thresh;
            if isfield(struct,'C')
                signifclus=struct.C<struct.thresh;
                t1=signifclus.*(1-signifpeak).*(127-struct.C/struct.thresh*126);
            else
                signifclus=0;
                t1=0;
            end
            t2=signifpeak.*(255-struct.P/struct.thresh*126);
            t3=(1-signifpeak).*(1-signifclus)*128;
            tt=(t1+t2+t3).*struct.mask*struct.thresh;

            trisurf(surf.tri,surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),...
                double(tt),'EdgeColor','none');
            set(gca,'Tag','SurfStatView 1 0');
            
            colorbar; colormap(spectral(256));
            cm=[zeros(1,3);
                zeros(127,1)   (0:126)'/127   ones(127,1); ones(1,3)*0.8;
                ones(127,1)    (0:126)'/127   zeros(127,1)];
            SurfStatColormap(cm);
            cb = SurfStatColLim( [0 255]*struct.thresh );

            set(cb,'YLim',[0 255]*struct.thresh);
            h=get(cb,'Children');
            set(h,'YData',[0 255]*struct.thresh);
            set(cb,'YTick',[1+(0:5)/5*126 129+(0:5)/5*126]*struct.thresh);
            pstr=strvcat(num2str(struct.thresh*(5:-1:1)'/5),'       0');
            set(cb,'YTickLabel',strvcat(pstr,pstr));
            xl=get(cb,'YLabel');
            set(xl,'String','P Cluster               P Vertex');
            
            dcm_obj=datacursormode(gcf);
            set(dcm_obj,'UpdateFcn',@SurfStatDataCursorQ,'DisplayStyle','window');
        end
        if isfield(struct,'Q')
            if ~isfield(struct,'thresh')
                struct.thresh=0.05;
            end

            t1=(struct.Q<struct.thresh).*(255-struct.Q/struct.thresh*253);
            t2=(struct.Q>=struct.thresh);
            tt=(t1+t2).*struct.mask*struct.thresh;

            trisurf(surf.tri,surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),...
                double(tt),'EdgeColor','none');
            set(gca,'Tag','SurfStatView 1 0');

            colorbar; colormap(spectral(256));
            cm=[zeros(1,3); ones(1,3)*0.8; ones(254,1) (0:253)'/254 zeros(254,1)];
            SurfStatColormap(cm);
            cb = SurfStatColLim( [0 255]*struct.thresh );

            set(cb,'YLim',[0 255]*struct.thresh);
            h=get(cb,'Children');
            set(h,'YData',[0 255]*struct.thresh);
            set(cb,'YTick',(2+(0:5)/5*253)*struct.thresh);
            set(cb,'YTickLabel',num2str(struct.thresh*(5:-1:0)'/5));
            xl=get(cb,'YLabel');
            set(xl,'String','Q');
            
            dcm_obj=datacursormode(gcf);
            set(dcm_obj,'UpdateFcn',@SurfStatDataCursorP,'DisplayStyle','window');
        end
    end
end

if isfield(surf,'lat')
    if ~isfield(surf,'vox')
        surf.vox=ones(1,3);
    end
    if ~isfield(surf,'origin');
        surf.origin=zeros(1,3);
    end    
    dim=size(surf.lat);
    x1=surf.origin(1)+((1:dim(1))-1)*surf.vox(1);
    y1=surf.origin(2)+((1:dim(2))-1)*surf.vox(2);
    z1=surf.origin(3)+((1:dim(3))-1)*surf.vox(3);
    [x,y,z]=meshgrid(x1,y1,z1);

    % Do a little bit of smoothing with FWHM=1.15 voxels;
    % since max(filt)>0.5, then mask>0.5 should still be the same as mask:
    [i,j,k]=ndgrid(-1:1,-1:1,-1:1);
    FWHM=1.15;
    filt=2.^(-4*(i.^2+j.^2+k.^2)/FWHM^2);
    filt=filt/sum(filt(:));

    if ~isstruct(struct)
        T=zeros(dim)+min(struct(:));
        T(surf.lat)=struct;
        T=permute(T,[2 1 3]);
        if (isempty(datathresh) & isempty(mask)) | ~isnan(xs) | ~isnan(ys) | ~isnan(zs) 
            if isnan(xs)
                xs=mean(x1);
            end
            if isnan(ys)
                ys=mean(y1);
            end
            if isnan(zs)
                zs=mean(z1);
            end
            s=slice(x,y,z,T,xs,ys,zs);
            set(s,'EdgeColor','none','FaceColor','interp');
            colorbar; colormap(spectral(256));
        end
        if ~isempty(datathresh)
            p1=patch(isosurface(x,y,z,T,datathresh));
            isonormals(x,y,z,T,p1);
            set(p1,'FaceColor',datacolor,'EdgeColor','none','FaceAlpha',dataalpha);
        end
    end

    if isfield(struct,'mask') | ~isempty(mask)
        M=zeros(dim);
        if isfield(struct,'mask') 
            M(surf.lat)=struct.mask;
        else
            M(surf.lat)=mask;
        end            
        M=convn(M,filt,'same');
        M=permute(M,[2 1 3]);
        p1=patch(isosurface(x,y,z,M,maskthresh));
        isonormals(x,y,z,M,p1);
        set(p1,'FaceColor',maskcolor,'EdgeColor','none','FaceAlpha',maskalpha);
    end

    if isfield(struct,'P')
        if ~isfield(struct,'thresh')
            struct.thresh=0.05;
        end
        P=ones(dim);
        P(surf.lat)=struct.P.*struct.mask+(1-struct.mask);
        P=permute(P,[2 1 3]);
        p2=patch(isosurface(x,y,z,P,struct.thresh));
        isonormals(x,y,z,P,p2)
        set(p2,'FaceColor',datacolor,'EdgeColor','none','FaceAlpha',dataalpha);

        if isfield(struct,'C')
            C=zeros(dim);
            C(surf.lat)=(struct.C<=struct.thresh).*struct.mask;
            C=convn(C,filt,'same');
            C=permute(C,[2 1 3]);
            p3=patch(isosurface(x,y,z,C,0.5));
            isonormals(x,y,z,C,p3)
            set(p3,'FaceColor',cluscolor,'EdgeColor','none','FaceAlpha',clusalpha);
        end
    end

    if isfield(struct,'Q')
        if ~isfield(struct,'thresh')
            struct.thresh=0.05;
        end
        Q=ones(dim);
        Q(surf.lat)=struct.Q.*struct.mask+(1-struct.mask);
        Q=permute(Q,[2 1 3]);
        p2=patch(isosurface(x,y,z,Q,struct.thresh));
        isonormals(x,y,z,Q,p2)
        set(p2,'FaceColor',datacolor,'EdgeColor','none','FaceAlpha',dataalpha);
    end
end

xlabel('x'); ylabel('y'); zlabel('z');
view(35,25);
daspect([1 1 1]); axis tight; axis vis3d; grid on;
lighting phong; material shiny; camlight; 

hold off;

title(title1);

whitebg(gcf,background);
set(gcf,'Color',background,'InvertHardcopy','off');

set(gcf,'PaperPosition',[0.25 2.5 6 4.5]);

return
end
