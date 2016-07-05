function qval = SurfStatQ( slm, mask );

%Q-values for False Discovey Rate of resels.
%
% Usage: qval = SurfStatQ( slm [, mask] );
%
% slm.t    = 1 x v vector of test statistics, v=#vertices.
% slm.df   = degrees of freedom.
% slm.dfs  = 1 x v vector of optional effective degrees of freedom.
% slm.k    = #variates.
% mask     = 1 x v logical vector, 1=inside, 0=outside, 
%          = ones(1,v), i.e. the whole surface, by default.
% The following are optional:
% slm.resl = e x k matrix of sum over observations of squares of
%           differences of normalized residuals along each edge.
% slm.tri  = t x 3 matrix of triangle indices, 1-based, t=#triangles.
% or
% slm.lat  = 3D logical array, 1=in, 0=out.
%
% qval.Q    = 1 x v vector of Q-values.
% qval.mask = copy of mask.

[l,v]=size(slm.t);
if nargin<2 | isempty(mask)
    mask=logical(ones(1,v));
end

df=zeros(2);
ndf=length(slm.df);
df(1,1:ndf)=slm.df;
df(2,1:2)=slm.df(ndf);
if isfield(slm,'dfs')
    df(1,ndf)=mean(slm.dfs(mask>0));
end

if isfield(slm,'du')
    [resels,reselspvert]=SurfStatResels(slm,mask);
else
    reselspvert=ones(1,v);
end
reselspvert=reselspvert(mask);

P_val=stat_threshold(0,1,0,df,[10 slm.t(1,mask)],[],[],[],slm.k,[],[],0);
P_val=P_val(2:length(P_val));
np=length(P_val);
[P_sort, index]=sort(P_val);
r_sort=reselspvert(index);
c_sort=cumsum(r_sort);
P_sort=P_sort./(c_sort+(c_sort<=0)).*(c_sort>0)*sum(r_sort);
m=1;
Q_sort=zeros(1,np);
for i=np:-1:1
    if P_sort(i)<m
        m=P_sort(i);
    end
    Q_sort(i)=m;
end
Q=zeros(1,np);
Q(index)=Q_sort;
qval.Q=ones(1,size(mask,2));
qval.Q(mask)=Q;
qval.mask=mask;

return
end

