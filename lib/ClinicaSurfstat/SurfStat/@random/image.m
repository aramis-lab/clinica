function image(m)
[sx,sv]=char(m);
[x,v]=double(m);
d1=0.66;
d2=0.88;

clf;
if ~isempty(x)
    [n,p]=size(x);
    minx=min(x);
    maxx=max(x);
    x=(x-ones(n,1)*minx)./(ones(n,1)*(maxx-minx+(maxx==minx))) ...
        + ones(n,1)*(maxx==minx);
    axes('position',[0.06 0.06 0.2 d1]);
    imagesc(x); colorbar; colormap(spectral(256));
    for k=1:p
        text(k,0,sx{k},'Rotation',90)
    end
    t=((0.97-0.06)/d1-1-1/2/n)*n;
    text(0.5,-t,'Mean');
end

if ~isempty(v)
    [n2,q]=size(v);
    n=sqrt(n2);
    q2=floor(sqrt(q-1))+1;
    t=0.06/(d2/q2-0.03)*n-0.5;
    s=0.015/(d2/q2-0.03)*n-0.5;
    for k=1:q
        [i,j]=ind2sub([q2,q2],k);
        axes('position',[0.28+(i-1)*d1/q2, 0.06+(q2-j)*d2/q2, ...
            d1/q2-0.03*3/4, d2/q2-0.03]);
        imagesc(reshape(v(:,k),n,n));
        axis off;
        text(n/2,-s,sv{k},'HorizontalAlignment','center')
        if k==1
            text(0.5,-t,'Variance');
        end
    end
end

background='white';
whitebg(gcf,background);
set(gcf,'Color',background,'InvertHardcopy','off');
set(gcf,'PaperPosition',[0.25 2.5 6 4.5]);
