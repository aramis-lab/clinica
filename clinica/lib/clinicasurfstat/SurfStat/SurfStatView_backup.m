function [ a, cb ] = SurfStatView( struct, surf, title, background);

%Viewer for surface data or P-value or Q-value structures.
% 
% Usage: [ a, cb ] = SurfStatView( struct, surf [,title [,background] ] );
% 
% struct        = 1 x v vector of data, v=#vertices, zeros(1,v) if empty,
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
% surf.coord    = 3 x v matrix of coordinates.
% surf.tri      = t x 3 matrix of triangle indices, 1-based, t=#triangles.
% title         = any string, data name by default.
% background    = background colour, any matlab ColorSpec, such as 
%   'white' (default), 'black'=='k', 'r'==[1 0 0], [1 0.4 0.6] (pink) etc.
%   Letter and line colours are inverted if background is dark (mean<0.5). 
%
% a  = vector of handles to the axes, left to right, top to bottom. 
% cb = handle to the colorbar.
%
% To change the colour limits, use SurfStatColLim( [min, max] ).
% To change the colour map, use e.g. SurfStatColormap( 'jet' ). 
% Surfaces can be edited in the figure window by clicking e.g. "Rotate 3D".
% If you want to customize the plot, modify the code in SurfStatViewData. 

if nargin<3 | isempty(title)
    title=inputname(1);
end
if nargin<4
    background='white'; % set up the default background color is white
end

if isempty(struct)
    struct=zeros(1,size(surf.coord,2));
end

if ~isstruct(struct) 
    [a,cb]=SurfStatViewData(struct,surf,title,background); % this is to view the data not in struct data type
else
    if isfield(struct,'P')  % This is to view the data from matlab struct type
        if ~isfield(struct,'thresh') % define the thresold for the map
            struct.thresh=0.05;
        end

        signifpeak=struct.P<struct.thresh; % Find the significant point in the map
        if isfield(struct,'C') % This is to define the corrected-pvalue in clinica
            signifclus=struct.C<struct.thresh;
            t1=signifclus.*(1-signifpeak).*(127-struct.C/struct.thresh*126);
            t2=signifpeak.*(255-struct.P/struct.thresh*126);
            t3=(1-signifpeak).*(1-signifclus)*128;
            tt=(t1+t2+t3).*struct.mask*struct.thresh;

            [a,cb]=SurfStatViewData(tt,surf,title,background);
            cm=[zeros(1,3);
                zeros(127,1)   (0:126)'/127   ones(127,1); ones(1,3)*0.8;
                ones(127,1)    (0:126)'/127   zeros(127,1)]; % This is the colormap that surfstat uses for prob value
            SurfStatColormap(cm); % the same with matlab colormap, this is just to regive the gcf a surfstat colormap
            cb = SurfStatColLim( [0 255]*struct.thresh );%  Sets the colour limits for SurfStatView. cb is the new colorbar

            set(cb,'XLim',[0 255]*struct.thresh);%Set object properties., here to set the x axis limitations
    % % %         h=get(cb,'Children'); % colorbar has no children, so this is an empty array
    % % %         set(h,'XData',[0 255]*struct.thresh); % see as below
            set(cb,'XTick',[1 64 127 129 192 255]*struct.thresh);% This is to set the xtick for the colorbar

            pstr1=num2str(round(struct.thresh*1000)/1000);
            pstr2=num2str(round(struct.thresh*1000/2)/1000);
            set(cb,'XTickLabel',strvcat(['     ' pstr1],['   ' pstr2],...
                ' ',['   0 ' pstr1],['   ' pstr2],'   0'));
            xl=get(cb,'XLabel');
            set(xl,'String','P Cluster               P Vertex');
            dcm_obj=datacursormode(gcf);
            set(dcm_obj,'UpdateFcn',@SurfStatDataCursorP,'DisplayStyle','window');
        else  % This is to define the uncorrected p value in clinica
            signifclus=0;
            t1=0;
            t2=signifpeak.*(255-struct.P/struct.thresh*126);
            t3=(1-signifpeak).*(1-signifclus)*128;
            tt=(t1+t2+t3).*struct.mask*struct.thresh;

            [a,cb]=SurfStatViewData(tt,surf,title,background);
            cm=[zeros(1,3);
                zeros(127,1)   (0:126)'/127   ones(127,1); ones(1,3)*0.8;
                ones(127,1)    (0:126)'/127   zeros(127,1)]; % This is the colormap that surfstat uses for prob value
            SurfStatColormap(cm); % the same with matlab colormap, this is just to regive the gcf a surfstat colormap
            cb = SurfStatColLim( [0 255]*struct.thresh );%  Sets the colour limits for SurfStatView. cb is the new colorbar

            set(cb,'XLim',[0 255]*struct.thresh);%Set object properties., here to set the x axis limitations
    % % %         h=get(cb,'Children'); % colorbar has no children, so this is an empty array
    % % %         set(h,'XData',[0 255]*struct.thresh); % see as below
            set(cb,'XTick',[1 64 127 129 192 255]*struct.thresh);% This is to set the xtick for the colorbar

            pstr1=num2str(round(struct.thresh*1000)/1000);
            pstr2=num2str(round(struct.thresh*1000/2)/1000);
            set(cb,'XTickLabel',strvcat(['     ' pstr1],['   ' pstr2],...
                ' ',['   0 ' pstr1],['   ' pstr2],'   0'));
            xl=get(cb,'XLabel');
            set(xl,'String','P Cluster               P Vertex');
            dcm_obj=datacursormode(gcf);
            set(dcm_obj,'UpdateFcn',@SurfStatDataCursorP,'DisplayStyle','window');
        end
    end
    if isfield(struct,'Q') % this is to define the F map in clinica
        if ~isfield(struct,'thresh')
            struct.thresh=0.05;
        end

        t1=(struct.Q<struct.thresh).*(255-struct.Q/struct.thresh*253);
        t2=(struct.Q>=struct.thresh);
        tt=(t1+t2).*struct.mask*struct.thresh;

        [a,cb]=SurfStatViewData(tt,surf,title,background);
        cm=[zeros(1,3); ones(1,3)*0.8; ones(254,1) (0:253)'/254 zeros(254,1)];
        SurfStatColormap(cm);
        cb = SurfStatColLim( [0 255]*struct.thresh );

        set(cb,'XLim',[0 255]*struct.thresh);
% % %         h=get(cb,'Children');
%         set(h,'XData',[0 255]*struct.thresh); % Xdata is just the
%         propertiy of an image, here cb is the colorbar, so maybe it is
%         cuz matlab advanced, before, XData can be used for colorbar
        set(cb,'XTick',(2+(0:5)/5*253)*struct.thresh);
        set(cb,'XTickLabel',num2str(struct.thresh*(5:-1:0)'/5));
        xl=get(cb,'XLabel');
        set(xl,'String','Q');
        dcm_obj=datacursormode(gcf);
        set(dcm_obj,'UpdateFcn',@SurfStatDataCursorQ,'DisplayStyle','window');
    end
end
    
return
end
