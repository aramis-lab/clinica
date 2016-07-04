function [ pcntvar, U, V ] = SurfStatPCA( Y, mask, X, k );

%Principal Components Analysis (PCA).
%
% Usage: [ pcntvar, U, V ] = SurfStatPCA( Y [,mask [,X [,k] ] ] );
%
% Y    = n x v matrix or n x v x k array of data, v=#vertices,
%        or memory map of same. 
% mask = 1 x v vector, 1=inside, 0=outside, default is ones(1,v),  
%        i.e. the whole surface.
% X    = model formula of type term, or scalar, or n x p design matrix of 
%        p covariates for the linear model. The PCA is done on the v x v 
%        correlations of the residuals and the components are standardized 
%        to have unit standard deviation about zero. If X=0, nothing is       
%        removed. If X=1, the mean (over rows) is removed (default).
% c    = number of components in PCA, default 4.
%
% pcntvar = 1 x c vector of percent variance explained by the components.
% U       = n x c matrix of components for the rows (observations).
% V       = c x v x k array of components for the columns (vertices).

maxchunk=2^20;
if isnumeric(Y)
    [n,v,k]=size(Y);
    isnum=true;
else
    Ym=Y;
    sz=Ym.Format{2};
    if length(sz)==2
        sz=sz([2 1]);
        k=1;
    else
        sz=sz([3 1 2]);
        k=sz(3);
    end
    n=sz(1);
    v=sz(2);
    isnum=false;
end    

if nargin<2 | isempty(mask)
    mask=ones(1,v);
end
if nargin<3 | isempty(X)
    X=1;
end
if nargin<4
    c=4;
end

if isa(X,'term')
    X=double(X);
elseif size(X,1)>1 & size(X,2)==1
    warning('If you don''t convert vectors to terms it can give unexpected results :-(') 
end
if size(X,1)==1
    X=repmat(X,n,1);
end
df=n-rank(X);

if isnum
    nc=1;
    chunk=v;
else
    nc=ceil(v*n*k/maxchunk);
    chunk=ceil(v/nc);
end
if ~isnum
    fprintf(1,'%s',[num2str(round(v*n*k*4/2^20)) ' Mb to PCA, % remaining: 100 ']);
    n10=floor(n/10);
end

A=zeros(n);
for ic=1:nc
    if ~isnum & rem(ic,n10)==0
        fprintf(1,'%s',[num2str(round(100*(1-ic/nc))) ' ']);
    end
    v1=1+(ic-1)*chunk;
    v2=min(v1+chunk-1,v);
    vc=v2-v1+1;
    if ~isnum
        if length(sz)==2
            Y=double(Ym.Data(1).Data(v1:v2,:)');
        else
            Y=double(permute(Ym.Data(1).Data(v1:v2,:,:),[3 1 2]));
        end
    end
    maskc=mask(v1:v2);
    if k==1
        Y=double(Y(:,maskc));
    else
        Y=double(reshape(Y(:,maskc,:),n,sum(maskc)*k));
    end
    if any(X(:)~=0)
        Y=Y-X*(pinv(X)*Y);
    end
    S=sum(Y.^2,1);
    Smhalf=(S>0)./sqrt(S+(S<=0));
    for i=1:n
        Y(i,:)=Y(i,:).*Smhalf;
    end
    A=A+Y*Y';
end
if ~isnum
    fprintf(1,'%s\n','Done');
    fprintf(1,'%s',[num2str(round(v*n*k*4/2^20)) ' Mb to PCA, % remaining: 100 ']);
end

[U, D]=eig(A);
[ds,is]=sort(-diag(D));
ds=-ds;
pcntvar=ds(1:c)'/sum(ds)*100;

U=U(:,is(1:c));
V=zeros(c,v*k);

if isnum
    V(:,repmat(mask,1,k))=U'*Y;
else
    for ic=1:nc
        if rem(ic,n10)==0
            fprintf(1,'%s',[num2str(round(100*(1-ic/nc))) ' ']);
        end
        v1=1+(ic-1)*chunk;
        v2=min(v1+chunk-1,v);
        vc=v2-v1+1;
        if length(sz)==2
            Y=double(Ym.Data(1).Data(v1:v2,:)');
        else
            Y=double(permute(Ym.Data(1).Data(v1:v2,:,:),[3 1 2]));
        end
        maskc=mask(v1:v2);
        if k==1
            Y=double(Y(:,maskc));
        else
            Y=double(reshape(Y(:,maskc,:),n,sum(maskc)*k));
        end
        if any(X(:)~=0)
            Y=Y-X*(pinv(X)*Y);
        end
        S=sum(Y.^2,1);
        Smhalf=(S>0)./sqrt(S+(S<=0));
        for i=1:n
            Y(i,:)=Y(i,:).*Smhalf;
        end
        maskcc=logical(zeros(1,v));
        maskcc(v1:v2)=maskc;
        if k>1
            maskcc=repmat(maskcc,1,k);
        end
        V(:,maskcc)=U'*Y;
    end
    fprintf(1,'%s\n','Done');
end

s=sign(abs(max(V,[],2))-abs(min(V,[],2)));
sv=sqrt(mean(V.^2,2));
V=diag(s./(sv+(sv<=0)).*(sv>0))*V;
U=U*diag(s*sqrt(df));
if k>1
    V=reshape(V,c,v,k);
end

return
end

