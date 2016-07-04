function [ pval, peak, clus, clusid ] = SurfStatP( slm, mask, clusthresh );

%Corrected P-values for vertices and clusters.
%
% Usage: [ pval, peak, clus, clusid ] = 
%             SurfStatP( slm [, mask, [, clusthresh] ] );
%
% slm.t       = l x v matrix of test statistics, v=#vertices; the first
%               row slm.t(1,:) is the test statistic, and the other 
%               rows are used to calculate cluster resels if k>1. See
%               SurfStatF for the precise definition of the extra rows.
% slm.df      = degrees of freedom.
% slm.dfs     = 1 x v vector of optional effective degrees of freedom.
% slm.k       = k=#variates.
% slm.resl    = e x k matrix of sum over observations of squares of
%               differences of normalized residuals along each edge.
% slm.tri     = 3 x t matrix of triangle indices, 1-based, t=#triangles.
% or
% slm.lat     = nx x ny x nz matrix, 1=in, 0=out, [nx,ny,nz]=size(volume). 
% mask        = 1 x v logical vector, 1=inside, 0=outside, v=#vertices, 
%             = ones(1,v), i.e. the whole surface, by default.
% clusthresh = P-value threshold or statistic threshold for 
%                   defining clusters, 0.001 by default.
%
% pval.P      = 1 x v vector of corrected P-values for vertices.
% pval.C      = 1 x v vector of corrected P-values for clusters.
% pval.mask   = copy of mask.
% peak.t      = np x 1 vector of peaks (local maxima).
% peak.vertid = np x 1 vector of vertex id's (1-based).
% peak.clusid = np x 1 vector of cluster id's that contain the peak.
% peak.P      = np x 1 vector of corrected P-values for the peak.
% clus.clusid = nc x 1 vector of cluster id numbers
% clus.nverts = nc x 1 vector of number of vertices in the cluster.
% clus.resels = nc x 1 vector of resels in the cluster.
% clus.P      = nc x 1 vector of corrected P-values for the cluster.
% clusid      =  1 x v vector of cluster id's for each vertex.
%
% Reference: Worsley, K.J., Andermann, M., Koulis, T., MacDonald, D. 
% & Evans, A.C. (1999). Detecting changes in nonisotropic images.
% Human Brain Mapping, 8:98-101.

[l,v]=size(slm.t);
if nargin<2 | isempty(mask)
    mask=logical(ones(1,v));
end
if nargin<3
    clusthresh=0.001;
end

df=zeros(2);
ndf=length(slm.df);
df(1,1:ndf)=slm.df;
df(2,1:2)=slm.df(ndf);
if isfield(slm,'dfs')
    df(1,ndf)=mean(slm.dfs(mask>0));
end

if v==1
    pval.P=stat_threshold(0,1,0,df,[10 slm.t(1)],[],[],[],slm.k,[],[],0);
    pval.P=pval.P(2);
    peak=[];
    clus=[];
    clusid=[];
    return
end

if clusthresh<1
    thresh=stat_threshold(0,1,0,df,clusthresh,[],[],[],slm.k,[],[],0);
else
    thresh=clusthresh;
end

[resels,reselspvert,edg]=SurfStatResels(slm,mask);

if max(slm.t(1,mask))<thresh
    pval.P=stat_threshold(resels,v,1,df,[10 slm.t],[],[],[],slm.k,[],[],0);
    pval.P=pval.P((1:v)+1);
    peak=[];
    clus=[];
    clusid=[];
else
    [peak,clus,clusid]=SurfStatPeakClus(slm,mask,thresh,reselspvert,edg);
    [pp,clpval]=stat_threshold(resels,v,1,df,...
        [10 peak.t' slm.t(1,:)],thresh,[10; clus.resels],[],slm.k,[],[],0);
    peak.P=pp((1:length(peak.t))+1)';
    pval.P=pp(length(peak.t)+(1:v)+1);
    if slm.k>1
        j=(slm.k-1):-2:0;
        sphere=zeros(1,slm.k);
        sphere(j+1)=exp((j+1)*log(2)+(j/2)*log(pi)+gammaln((slm.k+1)/2)- ...
            gammaln(j+1)-gammaln((slm.k+1-j)/2));
        sphere=sphere.*(4*log(2)).^(-(0:(slm.k-1))/2)/ndf;
        [pp,clpval]=stat_threshold(conv(resels,sphere),Inf,1,df,...
            [],thresh,[10; clus.resels],[],[],[],[],0);
    end
    clus.P=clpval(2:length(clpval));
    pval.C=interp1([0; clus.clusid],[1; clus.P],clusid);        
end
tlim=stat_threshold(resels,v,1,df,[0.5 1],[],[],[],slm.k,[],[],0);
tlim=tlim(2);
pval.P=pval.P.*(slm.t(1,:)>tlim)+(slm.t(1,:)<=tlim);
pval.mask=mask;

return
end

