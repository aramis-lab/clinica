function txt=SurfStatDataCursorP(empt,event_obj)
pos=get(event_obj,'Position');
h=get(event_obj,'Target');
v=get(h,'Vertices');
x=get(h,'FaceVertexCData');
id1=min(find(v(:,1)==pos(1)&v(:,2)==pos(2)&v(:,3)==pos(3)));
tag=get(get(h,'Parent'),'Tag');
[s,a,id0]=strread(tag,'%s %d %d');
id=id1+id0
c=get(get(h,'parent'),'CLim');
thresh=c(2)/255;
tt=x(id1)/thresh;
rtt=round(tt);
if rtt==128 | rtt==0
    p=NaN;
elseif rtt>=129
    p=(255-tt)/126*thresh;
else
    p=(127-tt)/126*thresh;
end
txt = {['x: ',num2str(pos(1))],...
       ['y: ',num2str(pos(2))],...
       ['z: ',num2str(pos(3))],...
       ['id: ',num2str(id)],...
       ['P: ',num2str(p)]};
return
end

