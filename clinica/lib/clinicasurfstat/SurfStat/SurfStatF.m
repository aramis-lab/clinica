function slm = SurfStatF( slm1, slm2 );

%F statistics for comparing two uni- or multi-variate fixed effects models.  
%
% Usage: slm = SurfStatF( slm1, slm2 );
%
% slm1 and slm2 must have the following fields:
% slm.X   = n x p design matrix.
% slm.df  = degrees of freedom.
% slm.SSE = k*(k+1)/2 x v matrix of sum of squares of errors.
%
% All fields of the bigger model are copied to slm, and 
% slm.k  = k=#variates.
% slm.df = [df1-df2, df2], where df1 and df2 are the min and max of slm1.df
%          and slm2.df, and SSE1 and SSE2 are the corresponding slm_.SSE's.
% slm.t  = l x v matrix of non-zero eigen values, in descending order, of 
%          F = (SSE1-SSE2)/(df1-df2)/(SSE2/df2), where l=min(k,df1-df2);  
%          slm.t(1,:) = Roy's maximum root = maximum F over all linear 
%          combinations of the k variates. k>3 is not programmed yet. 

if isfield(slm1,'r') | isfield(slm2,'r')
    warning('Mixed effects models not programmed yet.');
    return
end

if slm1.df>slm2.df
    X1=slm1.X;
    X2=slm2.X;
else
    X1=slm2.X;
    X2=slm1.X;
end

r=X1-X2*pinv(X2)*X1;
d=sum(r(:).^2)/sum(X1(:).^2);
if d>eps
    error('Models are not nested')
    return
end

if slm1.df>slm2.df
    df1=slm1.df(length(slm1.df));
    df2=slm2.df(length(slm2.df));
    SSE1=slm1.SSE;
    SSE2=slm2.SSE;
    slm=slm2;
else
    df1=slm2.df(length(slm2.df));
    df2=slm1.df(length(slm1.df));
    SSE1=slm2.SSE;
    SSE2=slm1.SSE;
    slm=slm1;
end

slm.df=[df1-df2, df2];
h=SSE1-SSE2;
if ndims(slm.coef)==2
    slm.k=1;
    slm.t=h./(SSE2+(SSE2<=0)).*(SSE2>0)*(df2/(df1-df2));
else
    [k2,v]=size(SSE2);
    k=round((sqrt(1+8*k2)-1)/2);
    slm.k=k;
    if k>3
        warning('Roy''s max root for k>3 not programmed yet.');
        return
    end
    l=min(k,df1-df2);
    slm.t=zeros(l,v);
    if k==2
        det=SSE2(1,:).*SSE2(3,:)-SSE2(2,:).^2;
        a11=SSE2(3,:).*h(1,:)-SSE2(2,:).*h(2,:);
        a21=SSE2(1,:).*h(2,:)-SSE2(2,:).*h(1,:);
        a12=SSE2(3,:).*h(2,:)-SSE2(2,:).*h(3,:);
        a22=SSE2(1,:).*h(3,:)-SSE2(2,:).*h(2,:);
        a0=a11.*a22-a12.*a21;
        a1=(a11+a22)/2;
        s1=real(sqrt(a1.^2-a0));
        d=(df2/(df1-df2))./(det+(det<=0)).*(det>0);
        slm.t(1,:)=(a1+s1).*d;
        if l==2
            slm.t(2,:)=(a1-s1).*d;
        end
    end
    if k==3
        det=SSE2(1,:).*(SSE2(3,:).*SSE2(6,:)-SSE2(5,:).^2) - ...
            SSE2(6,:).*SSE2(2,:).^2 + ...
            SSE2(4,:).*(SSE2(2,:).*SSE2(5,:)*2-SSE2(3,:).*SSE2(4,:));
        m1=SSE2(3,:).*SSE2(6,:)-SSE2(5,:).^2;
        m3=SSE2(1,:).*SSE2(6,:)-SSE2(4,:).^2;
        m6=SSE2(1,:).*SSE2(3,:)-SSE2(2,:).^2;
        m2=SSE2(4,:).*SSE2(5,:)-SSE2(2,:).*SSE2(6,:);
        m4=SSE2(2,:).*SSE2(5,:)-SSE2(3,:).*SSE2(4,:);
        m5=SSE2(2,:).*SSE2(4,:)-SSE2(1,:).*SSE2(5,:);
        a11=m1.*h(1,:)+m2.*h(2,:)+m4.*h(4,:);
        a12=m1.*h(2,:)+m2.*h(3,:)+m4.*h(5,:);
        a13=m1.*h(4,:)+m2.*h(5,:)+m4.*h(6,:);
        a21=m2.*h(1,:)+m3.*h(2,:)+m5.*h(4,:);
        a22=m2.*h(2,:)+m3.*h(3,:)+m5.*h(5,:);
        a23=m2.*h(4,:)+m3.*h(5,:)+m5.*h(6,:);
        a31=m4.*h(1,:)+m5.*h(2,:)+m6.*h(4,:);
        a32=m4.*h(2,:)+m5.*h(3,:)+m6.*h(5,:);
        a33=m4.*h(4,:)+m5.*h(5,:)+m6.*h(6,:);
        a0=-a11.*(a22.*a33-a23.*a32) ...
           +a12.*(a21.*a33-a23.*a31) ...
           -a13.*(a21.*a32-a22.*a31);
        a1=a22.*a33-a23.*a32 + ...
           a11.*a33-a13.*a31 + ...
           a11.*a22-a12.*a21;
        a2=-(a11+a22+a33);
        q=a1/3-a2.^2/9;
        r=(a1.*a2-3*a0)/6-a2.^3/27;
        s1=(r+sqrt(q.^3+r.^2)).^(1/3);
        z=zeros(3,v);
        z(1,:)=2*real(s1)-a2/3;
        z(2,:)=-real(s1)-a2/3+sqrt(3)*imag(s1);
        z(3,:)=-real(s1)-a2/3-sqrt(3)*imag(s1);
        z=sort(z,1,'descend');
        d=(df2/(df1-df2))./(det+(det<=0)).*(det>0);
        for j=1:l
            slm.t(j,:)=z(j,:).*d;
        end
    end
end

return
end

