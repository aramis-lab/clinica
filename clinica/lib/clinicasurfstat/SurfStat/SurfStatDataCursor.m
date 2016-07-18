function txt=SurfStatDataCursor(empt,event_obj)
pos=get(event_obj,'Position');
h=get(event_obj,'Target');
v=get(h,'Vertices');
x=get(h,'FaceVertexCData');
id1=min(find(v(:,1)==pos(1)&v(:,2)==pos(2)&v(:,3)==pos(3)));
tag=get(get(h,'Parent'),'Tag');
[s,a,id0]=strread(tag,'%s %d %d');
id=id1+id0
txt = {['x: ',num2str(pos(1))],...
       ['y: ',num2str(pos(2))],...
       ['z: ',num2str(pos(3))],...
       ['id: ',num2str(id)],...
       ['value: ',num2str(x(id1))]};
return
end

