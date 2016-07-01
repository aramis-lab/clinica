function [resels, reselspvert, edg] = SurfStatResels( slm, mask );

%Resels of surface or volume data inside a mask.
%
% Usage: [resels, reselspvert, edg] = SurfStatResels( slm [, mask] )
%
% slm.resl = e x k matrix of sum over observations of squares of
%            differences of normalized residuals along each edge.
% slm.tri  = t x 3 matrix of triangle indices, 1-based, t=#triangles.
% or
% slm.lat  = 3D logical array, 1=in, 0=out.
% mask     = 1 x v, 1=inside, 0=outside, v=#vertices,
%          = ones(1,v), i.e. the whole surface, by default.
%
% resels      = 1 x (D+1) vector of 0,...,D dimensional resels of the mask,
%             = EC of the mask if slm.resl is not given.
% reselspvert = 1 x v vector of D-dimensional resels per mask vertex.
% edg         = e x 2 matrix of edge indices, 1-based, e=#edges.

% The first version of this function (below) was much more elegant.
% There were two reasons why I had to abandon it: 
% 1) for volumetric data, the edg, tri and tet matrices become huge, often
% exceeding the memory capacity - they require 62 integers per voxel, as
% opposed to 12 integers per vertex for surface data. 
% 2) even if memory is not exceeded, the "unique" and "ismember" functions
% take a huge amount of time; the original version is ~7 times slower.
% To overcome 1), I had to resort to a slice-wise implementation for
% volumetric data. To cut down storage, the edg, tri and tet matrices are
% formed only for two adjacent slices in the volume. The slm.lat data is
% then copied into one slice, then the other, in an alternating fashion.
% Since the tetrahedral filling is also alternating, there is no need to
% recompute the edg, tri and tet matrices. However the ordering of the edg
% matrix is crucial - it must match the SurfStatEdg function so that 
% slm.resl is synchronised. The ordering of the tri and tet matrices is not
% important. 
% To avoid the "unique" function in 2), the edg and tri matrices were
% created by hand, rather than derived from the tet matrix. To avoid the
% "ismember" function, edge look-up was implemented by a sparse matrix. 
% This reduces execution time to ~18%. However it will not work for large
% dimensions, since the largest index of a sparse matrix is 2^31. If this
% is exceeded, then "interp1" is used which is slower but still reduces the
% time to ~45% that of "ismember". Traingle area look-up was avoided by
% simply re-computing trangle areas from the edges.
%
% Here is the original SurfStatResels function:
% 
% if isfield(surf,'tri')
%     tri=sort(surf.tri,2);
%     edg=unique([tri(:,[1 2]); tri(:,[1 3]); tri(:,[2 3])],'rows');
%     tet=[];
% end
% if isfield(surf,'lat')
%     % Fill lattice with 5 tetrahedra per cube
%     [I,J,K]=size(surf.lat);
%     [i,j,k]=ndgrid(1:int32(I),1:int32(J),1:int32(K));
%     i=i(:);
%     j=j(:);
%     k=k(:);
%     c1=find(rem(i+j+k,2)==0 & i<I & j<J & k<K);
%     c2=find(rem(i+j+k,2)==0 & i>1 & j<J & k<K);
%     clear i j k
%     IJ=I*J;
%     tet=[c1     c1+1    c1+1+I    c1+1+IJ;
%          c1     c1+I    c1+1+I    c1+I+IJ;
%          c1     c1+1+I  c1+1+IJ   c1+I+IJ;
%          c1     c1+IJ   c1+1+IJ   c1+I+IJ;
%          c1+1+I c1+1+IJ c1+I+IJ   c1+1+I+IJ;
%          c2-1   c2      c2-1+I    c2-1+IJ;
%          c2     c2-1+I  c2+I      c2+I+IJ;
%          c2     c2-1+I  c2-1+IJ   c2+I+IJ;
%          c2     c2-1+IJ c2+IJ     c2+I+IJ;
%          c2-1+I c2-1+IJ c2-1+I+IJ c2+I+IJ];
%     clear c1 c2
%     % Find triangles and edges
%     tri=unique([tet(:,[1 2 3]); tet(:,[1 2 4]); ...
%                 tet(:,[1 3 4]); tet(:,[2 3 4])],'rows');
%     edg=unique([tri(:,[1 2]); tri(:,[1 3]); tri(:,[2 3])],'rows');
%     % index by voxels in the lat
%     vid=cumsum(surf.lat(:)).*surf.lat(:);
%     % only inside the lat
%     edg=vid(edg(all(surf.lat(edg),2),:));
%     tri=vid(tri(all(surf.lat(tri),2),:));
%     tet=vid(tet(all(surf.lat(tet),2),:));
% end
%
% m=sum(mask(:));
% lkc(1,1)=m;
% %% LKC of edges
% maskedg=all(mask(edg),2);
% lkc(1,2)=sum(maskedg);
% if isfield(slm,'resl')
%     r1=mean(sqrt(slm.resl(maskedg,:)),2);
%     lkc(2,2)=sum(r1);
% end
% 
% %% LKC of triangles
% masktri=all(mask(tri),2);
% lkc(1,3)=sum(masktri);
% if isfield(slm,'resl')
%     [tf,loc]=ismember(tri(masktri,[1 2]),edg,'rows');  l12=slm.resl(loc,:);
%     [tf,loc]=ismember(tri(masktri,[1 3]),edg,'rows');  l13=slm.resl(loc,:);
%     [tf,loc]=ismember(tri(masktri,[2 3]),edg,'rows');  l23=slm.resl(loc,:);
%     a=max(4*l12.*l13-(l12+l13-l23).^2,0);
%     r2=mean(sqrt(a),2)/4;
%     lkc(2,3)=sum(mean(sqrt(l12)+sqrt(l13)+sqrt(l23),2))/2;
%     lkc(3,3)=sum(r2);
% end
% 
% if isfield(slm,'lat')
%     %% LKC of tetrahedra
%     masktet=all(mask(tet),2);
%     trimasktri=tri(masktri,:);
%     lkc(1,4)=sum(masktet);
%     if isfield(slm,'resl')
%         [tf,loc]=ismember(tet(masktet,[1 2 3]),trimasktri,'rows');  a4=a(loc,:);
%         [tf,loc]=ismember(tet(masktet,[1 2 4]),trimasktri,'rows');  a3=a(loc,:);
%         [tf,loc]=ismember(tet(masktet,[1 3 4]),trimasktri,'rows');  a2=a(loc,:);
%         [tf,loc]=ismember(tet(masktet,[2 3 4]),trimasktri,'rows');  a1=a(loc,:);
% 
%         [tf,loc]=ismember(tet(masktet,[1 2]),edg,'rows');  l12=slm.resl(loc,:);
%         [tf,loc]=ismember(tet(masktet,[1 3]),edg,'rows');  l13=slm.resl(loc,:);
%         [tf,loc]=ismember(tet(masktet,[2 3]),edg,'rows');  l23=slm.resl(loc,:);
%         [tf,loc]=ismember(tet(masktet,[1 4]),edg,'rows');  l14=slm.resl(loc,:);
%         [tf,loc]=ismember(tet(masktet,[2 4]),edg,'rows');  l24=slm.resl(loc,:);
%         [tf,loc]=ismember(tet(masktet,[3 4]),edg,'rows');  l34=slm.resl(loc,:);
% 
%         d12=4*l12.*l34-(l13+l24-l23-l14).^2;
%         d13=4*l13.*l24-(l12+l34-l23-l14).^2;
%         d14=4*l14.*l23-(l12+l34-l24-l13).^2;
% 
%         h=(a1<=0)|(a2<=0);
%         delta12=sum(mean(sqrt(l34).*pacos((d12-a1-a2)./sqrt(a1.*a2+h)/2.*(1-h)+h),2));
%         h=(a1<=0)|(a3<=0);
%         delta13=sum(mean(sqrt(l24).*pacos((d13-a1-a3)./sqrt(a1.*a3+h)/2.*(1-h)+h),2));
%         h=(a1<=0)|(a4<=0);
%         delta14=sum(mean(sqrt(l23).*pacos((d14-a1-a4)./sqrt(a1.*a4+h)/2.*(1-h)+h),2));
%         h=(a2<=0)|(a3<=0);
%         delta23=sum(mean(sqrt(l14).*pacos((d14-a2-a3)./sqrt(a2.*a3+h)/2.*(1-h)+h),2));
%         h=(a2<=0)|(a4<=0);
%         delta24=sum(mean(sqrt(l13).*pacos((d13-a2-a4)./sqrt(a2.*a4+h)/2.*(1-h)+h),2));
%         h=(a3<=0)|(a4<=0);
%         delta34=sum(mean(sqrt(l12).*pacos((d12-a3-a4)./sqrt(a3.*a4+h)/2.*(1-h)+h),2));
% 
%         r3=mean(sqrt(max((4*a1.*a2-(a1+a2-d12).^2)./(l34+(l34<=0)).*(l34>0),0)),2)/48;
% 
%         lkc(2,4)=(delta12+delta13+delta14+delta23+delta24+delta34)/(2*pi);
%         lkc(3,4)=sum(mean(sqrt(a1)+sqrt(a2)+sqrt(a3)+sqrt(a4),2))/8;
%         lkc(4,4)=sum(r3);
%     end
% end

if isfield(slm,'tri')
    tri=sort(slm.tri,2);
    edg=unique([tri(:,[1 2]); tri(:,[1 3]); tri(:,[2 3])],'rows');

    if nargin<2
        v=max(edg(:));
        mask=logical(zeros(1,v));
        mask(edg)=1;
    else
        if size(mask,1)>1
            mask=mask';
        end
        v=size(mask,2);
    end
        
    m=sum(mask(:));
    lkc(1,1)=m;
    %% LKC of edges
    maskedg=all(mask(edg),2);
    lkc(1,2)=sum(maskedg);
    if isfield(slm,'resl')
        r1=mean(sqrt(slm.resl(maskedg,:)),2);
        lkc(2,2)=sum(r1);
    end

    %% LKC of triangles
    masktri=all(mask(tri),2);
    lkc(1,3)=sum(masktri);
    if isfield(slm,'resl')
        [tf,loc]=ismember(tri(masktri,[1 2]),edg,'rows');  l12=slm.resl(loc,:);
        [tf,loc]=ismember(tri(masktri,[1 3]),edg,'rows');  l13=slm.resl(loc,:);
        [tf,loc]=ismember(tri(masktri,[2 3]),edg,'rows');  l23=slm.resl(loc,:);
        a=max(4*l12.*l13-(l12+l13-l23).^2,0);
        r2=mean(sqrt(a),2)/4;
        lkc(2,3)=sum(mean(sqrt(l12)+sqrt(l13)+sqrt(l23),2))/2;
        lkc(3,3)=sum(r2);
    end
    if nargout>=2
        reselspvert=zeros(v,1);
        for j=1:3
            reselspvert=reselspvert+accumarray(tri(masktri,j),r2,[v 1]);
        end
        D=2;
        reselspvert=reselspvert'/(D+1)/sqrt(4*log(2))^D;
    end
end
%%
if isfield(slm,'lat')
    edg=SurfStatEdg(slm);
    % The lattice is filled with 5 alternating tetrahedra per cube
    [I,J,K]=size(slm.lat);
    IJ=I*J;
    [i,j]=ndgrid(1:int32(I),1:int32(J));
    i=i(:);
    j=j(:);

    c1=int32(find(rem(i+j,2)==0 & i<I & j<J));
    c2=int32(find(rem(i+j,2)==0 & i>1 & j<J));
    c11=int32(find(rem(i+j,2)==0 & i==I & j<J));
    c21=int32(find(rem(i+j,2)==0 & i==I & j>1));
    c12=int32(find(rem(i+j,2)==0 & i<I & j==J));
    c22=int32(find(rem(i+j,2)==0 & i>1 & j==J));

    d1=find(rem(i+j,2)==1 & i<I & j<J)+IJ;
    d2=find(rem(i+j,2)==1 & i>1 & j<J)+IJ;

    tri1=[c1     c1+1    c1+1+I;  % bottom slice
        c1     c1+I    c1+1+I;
        c2-1   c2      c2-1+I;
        c2     c2-1+I  c2+I];
    tri2=[c1     c1+1    c1+1+IJ;  % between slices
        c1     c1+IJ   c1+1+IJ;
        c1     c1+I    c1+I+IJ;
        c1     c1+IJ   c1+I+IJ;
        c1     c1+1+I  c1+1+IJ;
        c1     c1+1+I  c1+I+IJ;
        c1     c1+1+IJ c1+I+IJ;
        c1+1+I c1+1+IJ c1+I+IJ;
        c2-1   c2      c2-1+IJ;
        c2     c2-1+IJ c2+IJ;
        c2-1   c2-1+I  c2-1+IJ;
        c2-1+I c2-1+IJ c2-1+I+IJ;
        c2     c2-1+I  c2+I+IJ;
        c2     c2-1+IJ c2+I+IJ;
        c2     c2-1+I  c2-1+IJ;
        c2-1+I c2-1+IJ c2+I+IJ;
        c11    c11+I    c11+I+IJ;
        c11    c11+IJ   c11+I+IJ;
        c21-I  c21      c21-I+IJ;
        c21    c21-I+IJ c21+IJ;
        c12    c12+1    c12+1+IJ;
        c12    c12+IJ   c12+1+IJ;
        c22-1  c22      c22-1+IJ;
        c22    c22-1+IJ c22+IJ];
    tri3=[d1     d1+1    d1+1+I;  % top slice
        d1     d1+I    d1+1+I;
        d2-1   d2      d2-1+I;
        d2     d2-1+I  d2+I];
    tet1=[c1     c1+1    c1+1+I    c1+1+IJ; % between slices
        c1     c1+I    c1+1+I    c1+I+IJ;
        c1     c1+1+I  c1+1+IJ   c1+I+IJ;
        c1     c1+IJ   c1+1+IJ   c1+I+IJ;
        c1+1+I c1+1+IJ c1+I+IJ   c1+1+I+IJ;
        c2-1   c2      c2-1+I    c2-1+IJ;
        c2     c2-1+I  c2+I      c2+I+IJ;
        c2     c2-1+I  c2-1+IJ   c2+I+IJ;
        c2     c2-1+IJ c2+IJ     c2+I+IJ;
        c2-1+I c2-1+IJ c2-1+I+IJ c2+I+IJ];

    v=sum(slm.lat(:));
    if nargin<2
        mask=logical(ones(1,v));
    else
        if size(mask,1)>1
            mask=mask';
        end
    end
    if nargout>=2
        reselspvert=zeros(v,1);
    end

    vs=cumsum(squeeze(sum(sum(slm.lat))));
    vs=int32([0; vs; vs(K)]);
    es=0;
    lat=logical(zeros(I,J,2));
    lat(:,:,1)=slm.lat(:,:,1);
    lkc=zeros(4);
    fprintf(1,'%s',[num2str(K) ' slices to resel, % remaining: 100 ']);
    n10=floor(K/10);
    for k=1:K
        if rem(k,n10)==0
            fprintf(1,'%s',[num2str(round(100*(1-k/K))) ' ']);
        end
        f=rem(k,2);
        if k<K
            lat(:,:,f+1)=slm.lat(:,:,k+1);
        else
            lat(:,:,f+1)=zeros(I,J);
        end
        vid=int32(cumsum(lat(:)).*lat(:))';
        if f
            edg1=edg(edg(:,1)>vs(k) & edg(:,1)<=vs(k+1),:)-vs(k);
            edg2=edg(edg(:,1)>vs(k) & edg(:,2)<=vs(k+2),:)-vs(k);
            tri=[vid(tri1(all(lat(tri1),2),:)); vid(tri2(all(lat(tri2),2),:))];
            mask1=mask((vs(k)+1):vs(k+2));
        else
            edg1=[edg(edg(:,1)>vs(k) & edg(:,2)<=vs(k+1),:)-vs(k)+vs(k+2)-vs(k+1);
                edg(edg(:,1)<=vs(k+1) & edg(:,2)>vs(k+1),2)-vs(k+1) ...
                edg(edg(:,1)<=vs(k+1) & edg(:,2)>vs(k+1),1)-vs(k)+vs(k+2)-vs(k+1)];
            edg2=[edg1; edg(edg(:,1)>vs(k+1) & edg(:,2)<=vs(k+2),:)-vs(k+1)];
            tri=[vid(tri3(all(lat(tri3),2),:)); vid(tri2(all(lat(tri2),2),:))];
            mask1=[mask((vs(k+1)+1):vs(k+2)) mask((vs(k)+1):vs(k+1))];
        end
        tet=vid(tet1(all(lat(tet1),2),:));

        m1=max(double(edg2(:,1)));
        ue=double(edg2(:,1))+m1*(double(edg2(:,2))-1);
        e=size(edg2,1);
        ae=1:e;
        if e<2^31
            sparsedg=sparse(ue,1,ae);
        end
        %%
        lkc1=zeros(4);
        lkc1(1,1)=sum(mask((vs(k)+1):vs(k+1)));
        
        %% LKC of edges
        maskedg=all(mask1(edg1),2);
        lkc1(1,2)=sum(maskedg);
        if isfield(slm,'resl')
            r1=mean(sqrt(slm.resl(find(maskedg)+es,:)),2);
            lkc1(2,2)=sum(r1);
        end
        
        %% LKC of triangles
        masktri=all(mask1(tri),2);
        lkc1(1,3)=sum(masktri);
        if isfield(slm,'resl')
            if e<2^31
                l12=slm.resl(sparsedg(tri(masktri,1)+m1*(tri(masktri,2)-1),1)+es,:);
                l13=slm.resl(sparsedg(tri(masktri,1)+m1*(tri(masktri,3)-1),1)+es,:);
                l23=slm.resl(sparsedg(tri(masktri,2)+m1*(tri(masktri,3)-1),1)+es,:);
            else
                l12=slm.resl(interp1(ue,ae,tri(masktri,1)+m1*(tri(masktri,2)-1),'nearest')+es,:);
                l13=slm.resl(interp1(ue,ae,tri(masktri,1)+m1*(tri(masktri,3)-1),'nearest')+es,:);
                l23=slm.resl(interp1(ue,ae,tri(masktri,2)+m1*(tri(masktri,3)-1),'nearest')+es,:);
            end
            a=max(4*l12.*l13-(l12+l13-l23).^2,0);
            r2=mean(sqrt(a),2)/4;
            lkc1(2,3)=sum(mean(sqrt(l12)+sqrt(l13)+sqrt(l23),2))/2;
            lkc1(3,3)=sum(r2);
            
            if nargout>=2 & K==1
                for j=1:3
                    if f
                        v1=tri(masktri,j)+vs(k);
                    else
                        v1=tri(masktri,j)+vs(k+1);
                        v1=v1-int32(v1>vs(k+2))*(vs(k+2)-vs(k));
                    end
                    reselspvert=reselspvert+accumarray(v1,r2,[v 1]);
                end
            end
        end
        
        %% LKC of tetrahedra
        masktet=all(mask1(tet),2);
        lkc1(1,4)=sum(masktet);
        if isfield(slm,'resl') & k<K
            if e<2^31
                l12=slm.resl(sparsedg(tet(masktet,1)+m1*(tet(masktet,2)-1),1)+es,:);
                l13=slm.resl(sparsedg(tet(masktet,1)+m1*(tet(masktet,3)-1),1)+es,:);
                l23=slm.resl(sparsedg(tet(masktet,2)+m1*(tet(masktet,3)-1),1)+es,:);
                l14=slm.resl(sparsedg(tet(masktet,1)+m1*(tet(masktet,4)-1),1)+es,:);
                l24=slm.resl(sparsedg(tet(masktet,2)+m1*(tet(masktet,4)-1),1)+es,:);
                l34=slm.resl(sparsedg(tet(masktet,3)+m1*(tet(masktet,4)-1),1)+es,:);
            else
                l12=slm.resl(interp1(ue,ae,tet(masktet,1)+m1*(tet(masktet,2)-1),'nearest')+es,:);
                l13=slm.resl(interp1(ue,ae,tet(masktet,1)+m1*(tet(masktet,3)-1),'nearest')+es,:);
                l23=slm.resl(interp1(ue,ae,tet(masktet,2)+m1*(tet(masktet,3)-1),'nearest')+es,:);
                l14=slm.resl(interp1(ue,ae,tet(masktet,1)+m1*(tet(masktet,4)-1),'nearest')+es,:);
                l24=slm.resl(interp1(ue,ae,tet(masktet,2)+m1*(tet(masktet,4)-1),'nearest')+es,:);
                l34=slm.resl(interp1(ue,ae,tet(masktet,3)+m1*(tet(masktet,4)-1),'nearest')+es,:);
            end
            a4=max(4*l12.*l13-(l12+l13-l23).^2,0);
            a3=max(4*l12.*l14-(l12+l14-l24).^2,0);
            a2=max(4*l13.*l14-(l13+l14-l34).^2,0);
            a1=max(4*l23.*l24-(l23+l24-l34).^2,0);

            d12=4*l12.*l34-(l13+l24-l23-l14).^2;
            d13=4*l13.*l24-(l12+l34-l23-l14).^2;
            d14=4*l14.*l23-(l12+l34-l24-l13).^2;

            h=(a1<=0)|(a2<=0);
            delta12=sum(mean(sqrt(l34).*pacos((d12-a1-a2)./sqrt(a1.*a2+h)/2.*(1-h)+h),2));
            h=(a1<=0)|(a3<=0);
            delta13=sum(mean(sqrt(l24).*pacos((d13-a1-a3)./sqrt(a1.*a3+h)/2.*(1-h)+h),2));
            h=(a1<=0)|(a4<=0);
            delta14=sum(mean(sqrt(l23).*pacos((d14-a1-a4)./sqrt(a1.*a4+h)/2.*(1-h)+h),2));
            h=(a2<=0)|(a3<=0);
            delta23=sum(mean(sqrt(l14).*pacos((d14-a2-a3)./sqrt(a2.*a3+h)/2.*(1-h)+h),2));
            h=(a2<=0)|(a4<=0);
            delta24=sum(mean(sqrt(l13).*pacos((d13-a2-a4)./sqrt(a2.*a4+h)/2.*(1-h)+h),2));
            h=(a3<=0)|(a4<=0);
            delta34=sum(mean(sqrt(l12).*pacos((d12-a3-a4)./sqrt(a3.*a4+h)/2.*(1-h)+h),2));

            r3=mean(sqrt(max((4*a1.*a2-(a1+a2-d12).^2)./(l34+(l34<=0)).*(l34>0),0)),2)/48;

            lkc1(2,4)=(delta12+delta13+delta14+delta23+delta24+delta34)/(2*pi);
            lkc1(3,4)=sum(mean(sqrt(a1)+sqrt(a2)+sqrt(a3)+sqrt(a4),2))/8;
            lkc1(4,4)=sum(r3);
            
            if nargout>=2
                for j=1:4
                    if f
                        v1=tet(masktet,j)+vs(k);
                    else
                        v1=tet(masktet,j)+vs(k+1);
                        v1=v1-int32(v1>vs(k+2))*(vs(k+2)-vs(k));
                    end
                    reselspvert=reselspvert+accumarray(v1,r3,[v 1]);
                end
            end
        end
        lkc=lkc+lkc1;
        es=es+size(edg1,1);
    end
    if nargout>=2
        D=2+(K>1);
        reselspvert=reselspvert'/(D+1)/sqrt(4*log(2))^D;
    end
    fprintf(1,'%s\n','Done');
end
%% resels
D=size(lkc,1)-1;
D2=size(lkc,2)-1;
tpltz=toeplitz((-1).^(0:D),(-1).^(0:D2));
lkcs=sum(tpltz.*lkc,2)';
lkcs=lkcs(1:max(find(abs(lkcs))));
resels=lkcs./sqrt(4*log(2)).^(0:D);

return
end

%%
function y=pacos(x)
y=acos(min(abs(x),1).*sign(x));
return
end
