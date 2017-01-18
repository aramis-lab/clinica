function [ a, cb ] = SurfStatViewData( data, surf, title, background );

%Basic viewer for surface data.
% 
% Usage: [ a, cb ] = SurfStatViewData( data, surf [,title [,background]] );
% 
% data        = 1 x v vector of data, v=#vertices
% surf.coord  = 3 x v matrix of coordinates.
% surf.tri    = t x 3 matrix of triangle indices, 1-based, t=#triangles.
% title       = any string, data name by default.
% background  = background colour, any matlab ColorSpec, such as 
%   'white' (default), 'black'=='k', 'r'==[1 0 0], [1 0.4 0.6] (pink) etc.
%   Letter and line colours are inverted if background is dark (mean<0.5).
%
% a  = vector of handles to the axes, left to right, top to bottom. 
% cb = handle to the colorbar.

if nargin<3 
    title=inputname(1);
end
if nargin<4
    background='white';
end

% find cut between hemispheres, assuming they are concatenated
t=size(surf.tri,1);
v=size(surf.coord,2);
tmax=max(surf.tri,[],2);
tmin=min(surf.tri,[],2);
% to save time, check that the cut is half way
if min(tmin(t/2+1:t))-max(tmax(1:t/2))==1
    cut=t/2;
    cuv=v/2;
else % check all cuts
    for i=1:t-1
        tmax(i+1)=max(tmax(i+1),tmax(i));
        tmin(t-i)=min(tmin(t-i),tmin(t-i+1));
    end
    cut=min([find((tmin(2:t)-tmax(1:t-1))==1) t]);
    cuv=tmax(cut);
end
tl=1:cut;
tr=(cut+1):t;
vl=1:cuv;
vr=(cuv+1):v;

clim=[min(data),max(data)];
if clim(1)==clim(2)
    clim=clim(1)+[-1 0]; % in case all the data is the same
end

% clf; clf deletes from the current figure all graphics objects whose 
% handles are not hidden (i.e., their HandleVisibility property is set to on).
% % % figure('Visible','off', 'position', [0, 0, 1650, 1050]);
figure('Visible','on', 'position', [0, 0, 1650, 1050]); % this is to
%display the image, make it visible
%here, we should find a way not to display the figure, but still got the
%same result that we want

colormap(spectral(256)); % this is the colormap that surfstat uses, set the
% figure's colormap to spectral(256)

h=0.39;
w=0.4;

r=max(surf.coord,[],2)-min(surf.coord,[],2);
w1=h/r(2)*r(1)*3/4;
h1=h/r(2)*r(1); % h/r(2)*r(3)

a(1)=axes('position',[0.055 0.62 h*3/4 w]);% defient gcf's axes propertities 
trisurf(surf.tri(tl,:),surf.coord(1,vl),surf.coord(2,vl),surf.coord(3,vl),...
    double(data(vl)),'EdgeColor','none'); % define the first trisurf image
view(-90,0); % This is to define the first trisurface displace angular 
daspect([1 1 1]); % sets the data aspect ratio.
axis tight; % Plot a surface. Set the axis limits to equal the range of the data so that the plot extends to the edges of the axes.
camlight; % define the direction where the light is from 
axis vis3d off; % not dispaly the axis
lighting phong; % control the lighting of SURFACE & PATCH objects, here, phong is a method for lighting 
material shiny; % sets the reflectance properties so that the object has a high specular reflectance relative to the diffuse and ambient light, and the color of the specular light depends only on the color of the light source
shading interp; % varies the color in each line segment and face by interpolating the colormap index or true color value across the line or face.

a(2)=axes('position',[0.3 0.58 w h]);
trisurf(surf.tri,surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),...
    double(data),'EdgeColor','none');
view(0,90); 
daspect([1 1 1]); axis tight; camlight; axis vis3d off;
lighting phong; material shiny; shading interp;

if cut<t
    a(3)=axes('position',[1-0.055-h*3/4 0.62 h*3/4 w]);
    trisurf(surf.tri(tr,:)-cuv,surf.coord(1,vr),surf.coord(2,vr),surf.coord(3,vr),...
        double(data(vr)),'EdgeColor','none');
    view(90,0);
    daspect([1 1 1]); axis tight; camlight; axis vis3d off;
    lighting phong; material shiny; shading interp;
else
    a(3)=axes('position',[1-0.055-h*3/4 0.62 h/r(2)*r(1)*3/4 w]);
    trisurf(surf.tri,surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),...
        double(data),'EdgeColor','none');
    view(180,0);
    daspect([1 1 1]); axis tight; camlight; axis vis3d off;
    lighting phong; material shiny; shading interp;
end

a(4)=axes('position',[0.055 0.29 h*3/4 w]);
trisurf(surf.tri(tl,:),surf.coord(1,vl),surf.coord(2,vl),surf.coord(3,vl),...
    double(data(vl)),'EdgeColor','none');
view(90,0); 
daspect([1 1 1]); axis tight; camlight; axis vis3d off;
lighting phong; material shiny; shading interp;

a(5)=axes('position',[0.3 0.18 w h]);
trisurf(surf.tri,surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),...
    double(data),'EdgeColor','none');
view(0,-90); 
daspect([1 1 1]); axis tight; camlight; axis vis3d off;
lighting phong; material shiny; shading interp;

if cut<t
    a(6)=axes('position',[1-0.055-h*3/4 0.29 h*3/4 w]);
    trisurf(surf.tri(tr,:)-cuv,surf.coord(1,vr),surf.coord(2,vr),surf.coord(3,vr),...
        double(data(vr)),'EdgeColor','none');
    view(-90,0);
    daspect([1 1 1]); axis tight; camlight; axis vis3d off;
    lighting phong; material shiny; shading interp;

    a(7)=axes('position',[0.055 0.02 w1 h1]);
    trisurf(surf.tri,surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),...
        double(data),'EdgeColor','none');
    view(180,0);
    daspect([1 1 1]); axis tight; camlight; axis vis3d off;
    lighting phong; material shiny; shading interp;

    a(8)=axes('position',[1-0.055-w1 0.03 w1 h1]);
    trisurf(surf.tri,surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),...
        double(data),'EdgeColor','none');
    view(0,0);
    daspect([1 1 1]); axis tight; camlight; axis vis3d off;
    lighting phong; material shiny; shading interp;
else
    a(6)=axes('position',[1-0.055-h*3/4 0.29 h/r(2)*r(1)*3/4 w]);
    trisurf(surf.tri,surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),...
        double(data),'EdgeColor','none');
    view(0,0);
    daspect([1 1 1]); axis tight; camlight; axis vis3d off;
    lighting phong; material shiny; shading interp;
end    
    
id0=[0 0 cuv 0 0 cuv 0 0];
for i=1:length(a)
    set(a(i),'CLim',clim);
    set(a(i),'Tag',['SurfStatView ' num2str(i) ' ' num2str(id0(i))]); %User-specified tag to associate with the axes, specified as a character vector. Tags provide a way to identify graphics objects. Use this property to find all objects with a specific tag within a plotting hierarchy, for example, searching for the tag using findobj.
end

cb=colorbar('location','South'); % define the colorbar location
set(cb,'Position',[0.35 0.085 0.3 0.03]); % define the colorbar specific position 
set(cb,'XAxisLocation','bottom'); % display the xaxis of the colorbar
% % % h=get(cb,'Title'); % get the Title property
% % % set(h,'String',title);

cb.Label.String = title;  % define the Label of the cb as title
cb.Label.FontSize = 12; % define the size of the fontsize of the title

whitebg(gcf,background);
set(gcf,'Color',background,'InvertHardcopy','off');

dcm_obj=datacursormode(gcf);
set(dcm_obj,'UpdateFcn',@SurfStatDataCursor,'DisplayStyle','window'); % This is to define the display style and updateFcn when you play with the 3D iamge when the figure 'Visible' is on.

set(gcf,'PaperPosition',[0.25 2.5 6 4.5]);

return
end
