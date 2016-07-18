function SurfStatColormap( map );

%Colormap function for SurfStatView.
%
% Usage: SurfStatColormap( map );
%
% Same as for matlab's colormap function - see help colormap.

colormap(map);
a=get(gcf,'Children');
k=0;
for i=1:length(a)
    tag=get(a(i),'Tag');
    if strcmp(tag,'Colorbar')
        cb=a(i);
    end
    if length(tag)>12 & strcmp(tag(1:12),'SurfStatView')
        k=k+1;
    end
end

if k>1
    set(cb,'location','South');
    set(cb,'Position',[0.35 0.085 0.3 0.03]);
    set(cb,'XAxisLocation','bottom');
end

return
end
