function image(t)
sx=char(t);
x=double(t);
if isempty(x)
    return
end
[n,p]=size(x);
minx=min(x);
maxx=max(x);
x=(x-ones(n,1)*minx)./(ones(n,1)*(maxx-minx+(maxx==minx))) ...
    + ones(n,1)*(maxx==minx);
d1=0.66;

clf;
axes('position',[0.06 0.06 0.2 d1]);
imagesc(x); colorbar; colormap(spectral(256)); 
for k=1:p
    text(k,0,sx{k},'Rotation',90)
end  
t=((0.97-0.06)/d1-1-1/2/n)*n;
text(0.5,-t,'Mean');

background='white';
whitebg(gcf,background);
set(gcf,'Color',background,'InvertHardcopy','off');
set(gcf,'PaperPosition',[0.25 2.5 6 4.5]);
