function slm = SurfStatT( slm, contrast );

%T statistics for a contrast in a univariate or multivariate model.
%
% Usage: slm = SurfStatT( slm, contrast );
%
% slm.X    = n x p design matrix.
% slm.V    = n x n x q variance matrix bases, normalised so that
%            mean(diag(slm.V))=1. If absent, assumes q=1 and slm.V=eye(n).
% slm.df   = degrees of freedom = n-rank(X).
% slm.coef = p x v x k matrix of coefficients of the linear model.
% slm.SSE  = k*(k+1)/2 x v matrix of sum of squares of errors.
% slm.r    = (q-1) x v matrix of coefficients of the first (q-1)
%            components of slm.V divided by their sum. 
% slm.dr   = (q-1) x 1 vector of increments in slm.r.
% contrast = n x 1 vector of contrasts in the observations, i.e.
%          = slm.X*slm.c', where slm.c is a contrast in slm.coef, or
%          = 1 x p vector of slm.c, padded with 0's if length(contrast)<p.
%
% slm.c   = 1 x p vector of contrasts in coefficents of the linear model.  
% slm.k   = k=#variates.
% slm.ef  = k x v matrix of effects.
% slm.sd  = k x v matrix of standard deviations of the effects.
% slm.t   = 1 x v vector of T = ef/sd if k=1, or Hotelling's T if k=2 or 3,
%           defined as the maximum T over all linear combinations of the k
%           variates; k>3 is not programmed yet.
% slm.dfs = 1 x v vector of effective degrees of freedom. Absent if q=1.
%
% Note that the contrast in the observations is used to determine the
% intended contrast in the model coefficients, slm.c. However there is some
% ambiguity in this when the model contains redundant terms. An example of
% such a model is 1 + Gender (Gender by itself does not contain redundant 
% terms). Only one of the ambiguous contrasts is estimable (i.e. has slm.sd
% < Inf), and this is the one chosen, though it may not be the contrast
% that you intended. To check this, compare the contrast in the 
% coefficients slm.c to the actual design matrix in slm.X. Note that the
% redundant columns of the design matrix have weights given by the rows of
% null(slm.X,'r')'

[n,p]=size(slm.X);
contrast=double(contrast);
pinvX=pinv(slm.X);
if length(contrast)<=p
    c=[contrast zeros(1,p-length(contrast))]';
    if sum((null(slm.X)'*c).^2)/sum(c.^2)>eps
        error('Contrast is not estimable :-(');
        return
    end
else
    c=pinvX*contrast;
    r=contrast-slm.X*c;
    if sum(r(:).^2)/sum(contrast(:).^2)>eps
%         warning('Contrast is not in the model :-(');
        error('Contrast is not in the model :-(');
    end
end
slm.c=c';

slm.df=slm.df(length(slm.df));
if ndims(slm.coef)==2
    slm.k=1;
    if ~isfield(slm,'r')
%% fixed effect
        if isfield(slm,'V')
            Vmh=inv(chol(slm.V)');
            pinvX=pinv(Vmh*slm.X);
        end
        Vc=sum((c'*pinvX).^2,2);
    else
%% mixed effect
        [q1,v]=size(slm.r);
        q=q1+1;
        
        nc=size(slm.dr,2);
        chunk=ceil(v/nc);
        irs=zeros(q1,v);
        for ic=1:nc
            v1=1+(ic-1)*chunk;
            v2=min(v1+chunk-1,v);
            vc=v2-v1+1;
            irs(:,v1:v2)=round(slm.r(:,v1:v2).*repmat(1./slm.dr(:,ic),1,vc));
        end
        [ur,ir,jr]=unique(irs','rows');
        nr=size(ur,1);
        slm.dfs=zeros(1,v);
        Vc=zeros(1,v);
        for ir=1:nr
            iv=(jr==ir);
            rv=mean(slm.r(:,iv),2);
            V=(1-sum(rv))*slm.V(:,:,q);
            for j=1:q1
                V=V+rv(j)*slm.V(:,:,j);
            end
            Vinv=inv(V);
            VinvX=Vinv*slm.X;
            Vbeta=pinv(slm.X'*VinvX);
            G=Vbeta*(VinvX');
            Gc=G'*c;
            R=Vinv-VinvX*G;
            E=zeros(q,1);
            for j=1:q
                E(j)=Gc'*slm.V(:,:,j)*Gc;
                RVV(:,:,j)=R*slm.V(:,:,j);
            end
            for j1=1:q
                for j2=j1:q
                    M(j1,j2)=sum(sum(RVV(:,:,j1).*(RVV(:,:,j2)')));
                    M(j2,j1)=M(j1,j2);
                end
            end
            vc=c'*Vbeta*c;
            iv=(ir==jr);
            Vc(iv)=vc;
            slm.dfs(iv)=vc^2/(E'*pinv(M)*E);
        end
    end
    slm.ef=c'*slm.coef;
    slm.sd=sqrt(Vc.*slm.SSE/slm.df);
    slm.t=slm.ef./(slm.sd+(slm.sd<=0)).*(slm.sd>0);
else
%% multivariate    
    [p,v,k]=size(slm.coef);
    slm.k=k;
    slm.ef=zeros(k,v);
    for j=1:k
        slm.ef(j,:)=c'*slm.coef(:,:,j);
    end
    j=1:k;
    jj=j.*(j+1)/2;
    vf=sum((c'*pinvX).^2,2)/slm.df;
    slm.sd=sqrt(vf*slm.SSE(jj,:));
    if k==2
        det=slm.SSE(1,:).*slm.SSE(3,:)-slm.SSE(2,:).^2;
        slm.t=slm.ef(1,:).^2.*slm.SSE(3,:) + ...
              slm.ef(2,:).^2.*slm.SSE(1,:) - ...
              2*slm.ef(1,:).*slm.ef(2,:).*slm.SSE(2,:);
    end
    if k==3
        det=slm.SSE(1,:).*(slm.SSE(3,:).*slm.SSE(6,:)-slm.SSE(5,:).^2) - ...
            slm.SSE(6,:).*slm.SSE(2,:).^2 + ...
            slm.SSE(4,:).*(slm.SSE(2,:).*slm.SSE(5,:)*2-slm.SSE(3,:).*slm.SSE(4,:));
        slm.t=      slm.ef(1,:).^2.*(slm.SSE(3,:).*slm.SSE(6,:)-slm.SSE(5,:).^2);
        slm.t=slm.t+slm.ef(2,:).^2.*(slm.SSE(1,:).*slm.SSE(6,:)-slm.SSE(4,:).^2);
        slm.t=slm.t+slm.ef(3,:).^2.*(slm.SSE(1,:).*slm.SSE(3,:)-slm.SSE(2,:).^2);
        slm.t=slm.t+2*slm.ef(1,:).*slm.ef(2,:).*(slm.SSE(4,:).*slm.SSE(5,:)-slm.SSE(2,:).*slm.SSE(6,:));
        slm.t=slm.t+2*slm.ef(1,:).*slm.ef(3,:).*(slm.SSE(2,:).*slm.SSE(5,:)-slm.SSE(3,:).*slm.SSE(4,:));
        slm.t=slm.t+2*slm.ef(2,:).*slm.ef(3,:).*(slm.SSE(2,:).*slm.SSE(4,:)-slm.SSE(1,:).*slm.SSE(5,:));
    end
    if k>3
         warning('Hotelling''s T for k>3 not programmed yet');
         return
    end
    slm.t=slm.t./(det+(det<=0)).*(det>0)/vf;
    slm.t=sqrt(slm.t+(slm.t<=0)).*(slm.t>0);
end

return
end



