function [ peak_threshold, extent_threshold, ...
    peak_threshold_1, extent_threshold_1, t, rho ] = ...
   stat_threshold( search_volume, num_voxels, fwhm, df, p_val_peak, ...
   cluster_threshold, p_val_extent, nconj, nvar, EC_file, expr, nprint );

%Thresholds and P-values of peaks and clusters of random fields in any D.
%
% [PEAK_THRESHOLD, EXTENT_THRESHOLD, PEAK_THRESHOLD_1 EXTENT_THRESHOLD_1] =
% STAT_THRESHOLD([ SEARCH_VOLUME [, NUM_VOXELS [, FWHM  [, DF [, P_VAL_PEAK 
% [, CLUSTER_THRESHOLD [, P_VAL_EXTENT [, NCONJ [, NVAR [, EC_FILE 
% [, TRANSFORM [, NPRINT ]]]]]]]]]]]] );
%
% If P-values are supplied, returns the thresholds for local maxima 
% or peaks (PEAK_THRESHOLD) (if thresholds are supplied then it
% returns P-values) for the following SPMs:
%  - Z, chi^2, t, F, Hotelling's T^2, Roy's maximum root, maximum 
%    canonical correlation, cross- and auto-correlation;
%  - scale space Z and chi^2; 
%  - conjunctions of NCONJ of these (except cross- and auto-correlation);
% and spatial extent of contiguous voxels above CLUSTER_THRESHOLD 
% (EXTENT_THRESHOLD) except conjunctions. Note that for the multivariate
% statistics spatial extent is measured in resels for the higher 
% dimensional space formed by adding sphere(s) to the search region.
%
% For peaks (local maxima), two methods are used:  
% random field theory (Worsley et al. 1996, Human Brain Mapping, 4:58-73),
% and a simple Bonferroni correction based on the number of voxels
% in the search region. The final threshold is the minimum of the two.
%
% For clusters, the method of Cao, 1999, Advances in Applied Probability,
% 31:577-593 is used. The cluster size is only accurate for large 
% CLUSTER_THRESHOLD (say >3) and large resels = SEARCH_VOLUME / FWHM^D. 
%
% PEAK_THRESHOLD_1 is the height of a single peak chosen in advance, and 
% EXTENT_THRESHOLD_1 is the extent of a single cluster chosen in advance
% of looking at the data, e.g. the nearest peak or cluster to a pre-chosen 
% voxel or ROI - see Friston KJ. Testing for anatomically specified 
% regional effects. Human Brain Mapping. 5(2):133-6, 1997.
%
% SEARCH_VOLUME is the volume of the search region in mm^3. The method for 
% finding PEAK_THRESHOLD works well for any value of SEARCH_VOLUME, even 
% SEARCH_VOLUME = 0 (default), which gives the threshold for the image at 
% a single voxel. The random field theory threshold is based on the 
% assumption that the search region is a sphere in 3D, which is a very 
% tight lower bound for any non-spherical region. 
% For a non-spherical D-dimensional region, set SEARCH_VOLUME to a vector  
% of the D+1 intrinsic volumes of the search region. In D=3 these are
% [Euler characteristic, 2 * caliper diameter, 1/2 surface area, volume].
% E.g. for a sphere of radius r in 3D, use [1, 4*r, 2*pi*r^2, 4/3*pi*r^3],
% which is equivalent to just SEARCH_VOLUME = 4/3*pi*r^3. 
% For a 2D search region, use [1, 1/2 perimeter length, area]. 
% For cross- and auto-correlations, where correlation is maximized over all
% pairs of voxels from two search regions, SEARCH_VOLUME has two rows
% interpreted as above - see Cao, J. & Worsley, K.J. (1999). The geometry 
% of correlation fields, with an application to functional connectivity of 
% the brain. Annals of Applied Probability, 9:1021-1057.
%
% NUM_VOXELS is the number of voxels (3D) or pixels (2D) in the search 
% volume. For cross- and auto-correlations, use two rows. Default is 1.
%
% FWHM is the fwhm in mm of a smoothing kernel applied to the data. Default
% is 0.0, i.e. no smoothing, which is roughly the case for raw fMRI data.
% For motion corrected fMRI data, use at least 6mm; for PET data, this would 
% be the scanner fwhm of about 6mm. Using the geometric mean of the three 
% x,y,z FWHM's is a very good approximation if they are different. 
% If you have already calculated resels, then set SEARCH_VOLUME=resels
% and FWHM=1. Cluster extents will then be in resels, rather than mm^3.
% For cross- and auto-correlations, use two rows. For scale space, use two 
% columns for the min and max fwhm (only for Z and chi^2 SPMs).  
%
% DF=[DF1 DF2; DFW1 DFW2 0] is a 2 x 2 matrix of degrees of freedom.
% If DF2 is 0, then DF1 is the df of the T statistic image.
% If DF1=Inf then it calculates thresholds for the Gaussian image. 
% If DF2>0 then DF1 and DF2 are the df's of the F statistic image.
% DFW1 and DFW2 are the numerator and denominator df of the FWHM image. 
% They are only used for calculating cluster resel p-values and thresholds
% (ignored if NVAR>1 or NCONJ>1 since there are no theoretical results).
% The random estimate of the local resels adds variability to the summed
% resels of the cluster. The distribution of the local resels summed over a 
% region depends on the unknown spatial correlation structure of the resels,
% so we assume that the resels are constant over the cluster, reasonable if the
% clusters are small. The cluster resels are then multiplied by the random resel 
% distribution at a point. This gives an upper bound to p-values and thresholds.
% If DF=[DF1 DF2] (row) then DFW1=DFW2=Inf, i.e. FWHM is fixed.
% If DF=[DF1; DFW1] (column) then DF2=0 and DFW2=DFW1, i.e. t statistic. 
% If DF=DF1 (scalar) then DF2=0 and DFW1=DFW2=Inf, i.e. t stat, fixed FWHM.
% If any component of DF >= 1000 then it is set to Inf. Default is Inf. 
%
% P_VAL_PEAK is the row vector of desired P-values for peaks. 
% If the first element is greater than 1, then they are
% treated as peak values and P-values are returned. Default is 0.05.
% To get P-values for peak values < 1, set the first element > 1,
% set the second to the deisired peak value, then discard the first result. 
%
% CLUSTER_THRESHOLD is the scalar threshold of the image for clusters.
% If it is <= 1, it is taken as a probability p, and the 
% cluter_thresh is chosen so that P( T > cluster_thresh ) = p. 
% Default is 0.001, i.e. P( T > cluster_thresh ) = 0.001. 
%
% P_VAL_EXTENT is the row vector of desired P-values for spatial extents
% of clusters of contiguous voxels above CLUSTER_THRESHOLD. 
% If the first element is greater than 1, then they are treated as 
% extent values and P-values are returned. Default is 0.05.
%
% NCONJ is the number of conjunctions. If NCONJ > 1, calculates P-values and 
% thresholds for peaks (but not clusters) of the minimum of NCONJ independent 
% SPM's - see Friston, K.J., Holmes, A.P., Price, C.J., Buchel, C.,
% Worsley, K.J. (1999). Multi-subject fMRI studies and conjunction analyses.
% NeuroImage, 10:385-396. Default is NCONJ = 1 (no conjunctions). 
%
% NVAR is the number of variables for multivariate equivalents of T and F 
% statistics, found by maximizing T^2 or F over all linear combinations of 
% variables, i.e. Hotelling's T^2 for DF1=1, Roy's maximum root for DF1>1. 
% For maximum canonical cross- and auto-correlations, use two rows. 
% See Worsley, K.J., Taylor, J.E., Tomaiuolo, F. & Lerch, J. (2004). Unified
% univariate and multivariate random field theory. Neuroimage, 23:S189-195.
% Default is 1, i.e. univariate statistics.
%
% EC_FILE: file for EC densities in IRIS Explorer latin module format.
% Ignored if empty (default).
%
% EXPR: matlab expression applied to threshold in t, e.g. 't=sqrt(t);'. 
%
% NPRINT: maximum number of P-values or thresholds to print. Default 5.
%
% Examples:
%
% T statistic, 1000000mm^3 search volume, 1000000 voxels (1mm voxel volume), 
% 20mm smoothing, 30 degrees of freedom, P=0.05 and 0.01 for peaks, cluster 
% threshold chosen so that P=0.001 at a single voxel, P=0.05 and 0.01 for extent:
%
%   stat_threshold(1000000,1000000,20,30,[0.05 0.01],0.001,[0.05 0.01]);
%
% peak_threshold = 5.0826    5.7853
% Cluster_threshold = 3.3852
% peak_threshold_1 = 4.8782    5.5955
% extent_threshold = 3340.3    6315.5
% extent_threshold_1 = 2911.5    5719.6
%
% Check: Suppose we are given peak values of 5.0826, 5.7853 and
% spatial extents of 3340.3, 6315.5 above a threshold of 3.3852,
% then we should get back our original P-values:
% 
%   stat_threshold(1000000,1000000,20,30,[5.0826 5.7853],3.3852,[3340.3 6315.5]);
%
% P_val_peak = 0.0500    0.0100
% P_val_peak_1 = 0.0318    0.0065
% P_val_extent = 0.0500    0.0100
% P_val_extent_1 = 0.0379    0.0074
%
% Another way of doing this is to use the fact that T^2 has an F 
% distribution with 1 and 30 degrees of freedom, and double the P values: 
% 
%   stat_threshold(1000000,1000000,20,[1 30],[0.1 0.02],0.002,[0.1 0.02]);
%
% peak_threshold = 25.8330   33.4702
% Cluster_threshold = 11.4596
% peak_threshold_1 = 20.7886   27.9741
% extent_threshold = 3297.8    6305.2
% extent_threshold_1 = 1936.4    4420.3
%
% Note that sqrt([25.8330 33.4702 11.4596])=[5.0826 5.7853 3.3852]
% in agreement with the T statistic thresholds, but the spatial extent
% thresholds are very close but not exactly the same.
%
% If the shape of the search region is known, then the 'intrinisic
% volume' must be supplied. E.g. for a spherical search region with
% a volume of 1000000 mm^3, radius r:
%
%   r=(1000000/(4/3*pi))^(1/3);
%   stat_threshold([1 4*r 2*pi*r^2 4/3*pi*r^3], ...
%                 1000000,20,30,[0.05 0.01],0.001,[0.05 0.01]);
% 
% gives the same results as our first example. A 100 x 100 x 100 cube 
%
%   stat_threshold([1 300 30000 1000000], ...
%                 1000000,20,30,[0.05 0.01],0.001,[0.05 0.01]);
% 
% gives slightly higher thresholds (5.0969, 5.7976) for peaks, but the cluster 
% thresholds are the same since they depend (so far) only on the volume and 
% not the shape. For a 2D circular search region of the same radius r, use
%
%   stat_threshold([1 pi*r pi*r^2],1000000,20,30,[0.05 0.01],0.001,[0.05 0.01]);
%
% A volume of 0 returns the usual uncorrected threshold for peaks,
% which can also be found by setting NUM_VOXELS=1, FWHM = 0:
%
%   stat_threshold(0,1,20,30,[0.05 0.01])
%   stat_threshold(1,1, 0,30,[0.05 0.01]);
%
% For non-isotropic fields, replace the SEARCH_VOLUME by the vector of resels 
% (see Worsley et al. 1996, Human Brain Mapping, 4:58-73), and set FWHM=1, 
% so that EXTENT_THRESHOLD's are measured in resels, not mm^3. If the 
% resels are estimated from residuals, add an extra row to DF (see above).
%
% For the conjunction of 2 and 3 T fields as above:
%
%   stat_threshold(1000000,1000000,20,30,[0.05 0.01],[],[],2)
%   stat_threshold(1000000,1000000,20,30,[0.05 0.01],[],[],3)
%
% returns lower thresholds of [3.0251 3.3984] and [2.2336 2.5141] resp. Note
% that there are as yet no theoretical results for extents so they are NaN.
% 
% For Hotelling's T^2 e.g. for linear models of deformations (NVAR=3):
%
%   stat_threshold(1000000,1000000,20,[1 30],[0.05 0.01],[],[],1,3)
%
% For Roy's max root, e.g. for effective connectivity using deformations,
%
%   stat_threshold(1000000,1000000,20,[3 30],[0.05 0.01],[],[],1,3)
%
% There are no theoretical results for mulitvariate extents so they are NaN.
% Check: A Hotelling's T^2 with Inf df is the same as a chi^2 with NVAR df 
%
%   stat_threshold(1000000,1000000,20,[1 Inf],[],[],[],[],NVAR)
%   stat_threshold(1000000,1000000,20,[NVAR Inf])*NVAR
%
% should give the same answer!
%
% Cross- and auto-correlations: to threshold the auto-correlation of fmrilm
% residuals (100 df) over all pairs of 1000000 voxels in a 1000000mm^3 
% region with 20mm FWHM smoothing, use df=99 (one less for the correlation): 
%
%   stat_threshold([1000000; 1000000], [1000000; 1000000], 20, 99, 0.05)
%
% which gives a T statistic threshold of 6.7923 with 99 df. To convert to 
% a correlation, use 6.7923/sqrt(99+6.7923^2)=0.5638. Note that auto-
% correlations are symmetric, so the P-value is in fact 0.05/2=0.025; on the
% other hand, we usually want a two sided test of both positive and negative
% correlations, so the P-value should be doubled back to 0.05. For cross-
% correlations, e.g. between n=321 GM density volumes (1000000mm^3 ball, 
% FWHM=13mm) and cortical thickness on the average surface (250000mm^2, 
% FWHM=18mm), use 321-2=319 df for removing the constant and the correlation.
% Note that we use P=0.025 to threshold both +ve and -ve correlations.
%
%   r=(1000000/(4/3*pi))^(1/3);
%   stat_threshold([1 4*r 2*pi*r^2 4/3*pi*r^3; 1 0 250000 0], ...
%      [1000000; 40962], [13; 18], 319, 0.025)
%
% This gives a T statistic threshold of 6.7158 with 319 df. To convert to 
% a correlation, use 6.7158/sqrt(319+6.7158^2)=0.3520. Note that NVAR can be
% greater than one for maximum canonical cross- and auto-correlations 
% between all pairs of multivariate imaging data. Finally, the spatial
% extent for 5D clusters of connected correlations is 5.9228e+006 mm^5. 
%
% Scale space: suppose we smooth 1000000 voxels in a 1000000mm^3 region
% from 6mm to 30mm in 10 steps ( e.g. fwhm=6*(30/6).^((0:9)/9) ):
%
%   stat_threshold(1000000, 1000000*10, [6 30])
%
% gives a scale space threshold of 5.0663, only slightly higher than the 
% 5.0038 you get from using the highest resolution data only, i.e.
%
%   stat_threshold(1000000, 1000000, 6)

%############################################################################
% COPYRIGHT:   Copyright 2003 K.J. Worsley 
%              Department of Mathematics and Statistics,
%              McConnell Brain Imaging Center, 
%              Montreal Neurological Institute,
%              McGill University, Montreal, Quebec, Canada. 
%              keith.worsley@mcgill.ca , www.math.mcgill.ca/keith
%
%              Permission to use, copy, modify, and distribute this
%              software and its documentation for any purpose and without
%              fee is hereby granted, provided that this copyright
%              notice appears in all copies. The author and McGill University
%              make no representations about the suitability of this
%              software for any purpose.  It is provided "as is" without
%              express or implied warranty.
%############################################################################

% Defaults:
if nargin<1;  search_volume=[];  end
if nargin<2;  num_voxels=[];  end
if nargin<3;  fwhm=[];  end
if nargin<4;  df=[];  end
if nargin<5;  p_val_peak=[];  end
if nargin<6;  cluster_threshold=[];  end
if nargin<7;  p_val_extent=[];  end
if nargin<8;  nconj=[];  end
if nargin<9;  nvar=[];  end
if nargin<10;  EC_file=[];  end
if nargin<11;  expr=[];  end
if nargin<12;  nprint=[];  end

if isempty(search_volume);  search_volume=0;  end
if isempty(num_voxels);  num_voxels=1;  end
if isempty(fwhm);  fwhm=0.0;  end
if isempty(df);  df=Inf;  end
if isempty(p_val_peak);  p_val_peak=0.05;  end
if isempty(cluster_threshold);  cluster_threshold=0.001;  end
if isempty(p_val_extent);  p_val_extent=0.05;  end
if isempty(nconj);  nconj=1;  end
if isempty(nvar);  nvar=1;  end
if isempty(nprint);  nprint=5;  end

if size(fwhm,1)==1; fwhm(2,:)=fwhm; end
if size(fwhm,2)==1; scale=1; else scale=fwhm(1,2)/fwhm(1,1); fwhm=fwhm(:,1); end;
isscale=(scale>1); 

if length(num_voxels)==1; num_voxels(2,1)=1; end

if size(search_volume,2)==1
   radius=(search_volume/(4/3*pi)).^(1/3);
   search_volume=[ones(length(radius),1) 4*radius 2*pi*radius.^2 search_volume];
end
if size(search_volume,1)==1
   search_volume=[search_volume; [1 zeros(1,size(search_volume,2)-1)]];
end
lsv=size(search_volume,2);
fwhm_inv=all(fwhm>0)./(fwhm+any(fwhm<=0));
resels=search_volume.*repmat(fwhm_inv,1,lsv).^repmat(0:lsv-1,2,1);
invol=resels.*(4*log(2)).^(repmat(0:lsv-1,2,1)/2);
for k=1:2
   D(k,1)=max(find(invol(k,:)))-1;
end

% determines which method was used to estimate fwhm (see fmrilm or multistat): 
df_limit=4;

if length(df)==1; df=[df 0]; end
if size(df,1)==1; df=[df; Inf Inf; Inf Inf]; end
if size(df,2)==1; df=[df df]; df(1,2)=0; end
if size(df,1)==2; df=[df; df(2,:)]; end

% is_tstat=1 if it is a t statistic
is_tstat=(df(1,2)==0);
if is_tstat
   df1=1;
   df2=df(1,1);
else
   df1=df(1,1);
   df2=df(1,2);
end
if df2 >= 1000; df2=Inf; end
df0=df1+df2;

dfw1=df(2:3,1);
dfw2=df(2:3,2);
for k=1:2
    if dfw1(k) >= 1000; dfw1(k)=Inf; end
    if dfw2(k) >= 1000; dfw2(k)=Inf; end
end

if length(nvar)==1; nvar(2,1)=df1; end

if isscale & (D(2)>1 | nvar(1,1)>1 | df2<Inf)
   D
   nvar
   df2
   fprintf('Cannot do scale space.');
   return
end

Dlim=D+[scale>1; 0];
DD=Dlim+nvar-1;

% Values of the F statistic:
t=((1000:-1:1)'/100).^4;
% Find the upper tail probs cumulating the F density using Simpson's rule:
if df2==Inf
   u=df1*t;
   b=exp(-u/2-log(2*pi)/2+log(u)/4)*df1^(1/4)*4/100;
else  
   u=df1*t/df2;
   b=exp(-df0/2*log(1+u)+log(u)/4-betaln(1/2,(df0-1)/2))*(df1/df2)^(1/4)*4/100;
end
t=[t; 0];
b=[b; 0];
n=length(t);
sb=cumsum(b);
sb1=cumsum(b.*(-1).^(1:n)');
pt1=sb+sb1/3-b/3;
pt2=sb-sb1/3-b/3;
tau=zeros(n,DD(1)+1,DD(2)+1);
tau(1:2:n,1,1)=pt1(1:2:n);
tau(2:2:n,1,1)=pt2(2:2:n);
tau(n,1,1)=1;
tau(:,1,1)=min(tau(:,1,1),1);

% Find the EC densities:
u=df1*t;
for d=1:max(DD)
   for e=0:min(min(DD),d)
      s1=0;
      cons=-((d+e)/2+1)*log(pi)+gammaln(d)+gammaln(e+1);
      for k=0:(d-1+e)/2
         [i,j]=ndgrid(0:k,0:k);
         if df2==Inf
            q1=log(pi)/2-((d+e-1)/2+i+j)*log(2);
         else
            q1=(df0-1-d-e)*log(2)+gammaln((df0-d)/2+i)+gammaln((df0-e)/2+j) ...
               -gammalni(df0-d-e+i+j+k)-((d+e-1)/2-k)*log(df2);
         end
         q2=cons-gammalni(i+1)-gammalni(j+1)-gammalni(k-i-j+1) ...
            -gammalni(d-k-i+j)-gammalni(e-k-j+i+1);
         s2=sum(sum(exp(q1+q2)));
         if s2>0
            s1=s1+(-1)^k*u.^((d+e-1)/2-k)*s2;
         end
      end
      if df2==Inf
         s1=s1.*exp(-u/2);
      else
         s1=s1.*exp(-(df0-2)/2*log(1+u/df2));
      end
      if DD(1)>=DD(2)
         tau(:,d+1,e+1)=s1;
         if d<=min(DD)
            tau(:,e+1,d+1)=s1;
         end
      else
         tau(:,e+1,d+1)=s1;      
         if d<=min(DD)
            tau(:,d+1,e+1)=s1;
         end
      end
   end
end

% For multivariate statistics, add a sphere to the search region:
a=zeros(2,max(nvar));
for k=1:2
   j=(nvar(k)-1):-2:0;
   a(k,j+1)=exp(j*log(2)+j/2*log(pi) ...
      +gammaln((nvar(k)+1)/2)-gammaln((nvar(k)+1-j)/2)-gammaln(j+1));
end
rho=zeros(n,Dlim(1)+1,Dlim(2)+1);
for k=1:nvar(1)
   for l=1:nvar(2)
      rho=rho+a(1,k)*a(2,l)*tau(:,(0:Dlim(1))+k,(0:Dlim(2))+l);
   end
end

if is_tstat
    if nvar==1
        t=[sqrt(t(1:(n-1))); -flipdim(sqrt(t),1)];
        rho=[rho(1:(n-1),:,:); flipdim(rho,1)]/2;
        for i=0:D(1)
            for j=0:D(2)
                rho(n-1+(1:n),i+1,j+1)=-(-1)^(i+j)*rho(n-1+(1:n),i+1,j+1);
            end
        end
        rho(n-1+(1:n),1,1)=rho(n-1+(1:n),1,1)+1;
        n=2*n-1;
    else
        t=sqrt(t);
    end
end

% For scale space:
if scale>1
   kappa=D(1)/2;
   tau=zeros(n,D(1)+1);
   for d=0:D(1)
      s1=0;
      for k=0:d/2
         s1=s1+(-1)^k/(1-2*k)*exp(gammaln(d+1)-gammaln(k+1)-gammaln(d-2*k+1) ...
            +(1/2-k)*log(kappa)-k*log(4*pi))*rho(:,d+2-2*k,1);
      end
      if d==0
         cons=log(scale);
      else
         cons=(1-1/scale^d)/d;
      end
      tau(:,d+1)=rho(:,d+1,1)*(1+1/scale^d)/2+s1*cons;
   end
   rho(:,1:(D(1)+1),1)=tau;
end

if D(2)==0
   d=D(1);
   if nconj>1
      % Conjunctions:
      b=gamma(((0:d)+1)/2)/gamma(1/2);
      for i=1:d+1
         rho(:,i,1)=rho(:,i,1)/b(i);
      end
      m1=zeros(n,d+1,d+1);
      for i=1:d+1
         j=i:d+1;
         m1(:,i,j)=rho(:,j-i+1,1);
      end
      for k=2:nconj
         for i=1:d+1
            for j=1:d+1
               m2(:,i,j)=sum(rho(:,1:d+2-i,1).*m1(:,i:d+1,j),2);
            end
         end
         m1=m2;
      end
      for i=1:d+1
         rho(:,i,1)=m1(:,1,i)*b(i);
      end
   end
   
   if ~isempty(EC_file)
      if d<3
         rho(:,(d+2):4,1)=zeros(n,4-d-2+1);
      end
      fid=fopen(EC_file,'w');
      % first 3 are dimension sizes as 4-byte integers:
      fwrite(fid,[n max(d+2,5) 1],'int');
      % next 6 are bounding box as 4-byte floats: 
      fwrite(fid,[0 0 0; 1 1 1],'float');
      % rest are the data as 4-byte floats:
      if ~isempty(expr)
         eval(expr);
      end
      fwrite(fid,t,'float');
      fwrite(fid,rho,'float');
      fclose(fid);
   end
end

if all(fwhm>0)
   pval_rf=zeros(n,1);
   for i=1:D(1)+1
      for j=1:D(2)+1
         pval_rf=pval_rf+invol(1,i)*invol(2,j)*rho(:,i,j);
      end
   end
else
   pval_rf=Inf;
end

% Bonferroni 
pt=rho(:,1,1);
pval_bon=abs(prod(num_voxels))*pt;

% Minimum of the two:
pval=min(pval_rf,pval_bon);

tlim=1;
if p_val_peak(1) <= tlim
   peak_threshold=minterp1(pval,t,p_val_peak);
   if length(p_val_peak)<=nprint
      peak_threshold
   end
else
   % p_val_peak is treated as a peak value:
   P_val_peak=interp1(t,pval,p_val_peak);
   i=isnan(P_val_peak);
   P_val_peak(i)=(is_tstat & (p_val_peak(i)<0));
   peak_threshold=P_val_peak;
   if length(p_val_peak)<=nprint
      P_val_peak
   end
end

if all(fwhm<=0) | any(num_voxels<0)
   peak_threshold_1=p_val_peak+NaN;
   extent_threshold=p_val_extent+NaN;
   extent_threshold_1=extent_threshold;
   return
end

% Cluster_threshold:

if cluster_threshold > tlim
   tt=cluster_threshold;
else
   % cluster_threshold is treated as a probability:
   tt=minterp1(pt,t,cluster_threshold);
   if nprint>0 
       Cluster_threshold=tt
   end
end

d=sum(D);
rhoD=interp1(t,rho(:,D(1)+1,D(2)+1),tt);
p=interp1(t,pt,tt);

% Pre-selected peak:

pval=rho(:,D(1)+1,D(2)+1)./rhoD;

if p_val_peak(1) <= tlim 
   peak_threshold_1=minterp1(pval,t, p_val_peak);
   if length(p_val_peak)<=nprint
      peak_threshold_1
   end
else
   % p_val_peak is treated as a peak value:
   P_val_peak_1=interp1(t,pval,p_val_peak);
   i=isnan(P_val_peak_1);
   P_val_peak_1(i)=(is_tstat & (p_val_peak(i)<0));
   peak_threshold_1=P_val_peak_1;
   if length(p_val_peak)<=nprint
      P_val_peak_1
   end
end

if  d==0 | nconj>1 | nvar(1)>1 | scale>1
    extent_threshold=p_val_extent+NaN;
    extent_threshold_1=extent_threshold;
    if length(p_val_extent)<=nprint
       extent_threshold
       extent_threshold_1
    end
    return
end

% Expected number of clusters:

EL=invol(1,D(1)+1)*invol(2,D(2)+1)*rhoD;
cons=gamma(d/2+1)*(4*log(2))^(d/2)/fwhm(1)^D(1)/fwhm(2)^D(2)*rhoD/p;

if df2==Inf & dfw1(1)==Inf & dfw1(2)==Inf
   if p_val_extent(1) <= tlim 
      pS=-log(1-p_val_extent)/EL;
      extent_threshold=(-log(pS)).^(d/2)/cons;
      pS=-log(1-p_val_extent);
      extent_threshold_1=(-log(pS)).^(d/2)/cons;
      if length(p_val_extent)<=nprint
         extent_threshold
         extent_threshold_1
      end
   else
      % p_val_extent is now treated as a spatial extent:
      pS=exp(-(p_val_extent*cons).^(2/d));
      P_val_extent=1-exp(-pS*EL);
      extent_threshold=P_val_extent;
      P_val_extent_1=1-exp(-pS);
      extent_threshold_1=P_val_extent_1;
      if length(p_val_extent)<=nprint
         P_val_extent
         P_val_extent_1
      end
   end
else
   % Find dbn of S by taking logs then using fft for convolution:
   ny=2^12;
   a=d/2;
   b2=a*10*max(sqrt(2/(min(df1+df2,min(dfw1)))),1);
   if df2<Inf
      b1=a*log((1-(1-0.000001)^(2/(df2-d)))*df2/2);
   else
      b1=a*log(-log(1-0.000001));
   end
   dy=(b2-b1)/ny;
   b1=round(b1/dy)*dy;
   y=((1:ny)'-1)*dy+b1;
   numrv=1+(d+(D(1)>0)+(D(2)>0))*(df2<Inf)+...
       (D(1)*(dfw1(1)<Inf)+(dfw2(1)<Inf))*(D(1)>0)+...;
       (D(2)*(dfw1(2)<Inf)+(dfw2(2)<Inf))*(D(2)>0);
   f=zeros(ny,numrv);
   mu=zeros(1,numrv);
   if df2<Inf
      % Density of log(Beta(1,(df2-d)/2)^(d/2)):
      yy=exp(y./a)/df2*2;  
      yy=yy.*(yy<1);
      f(:,1)=(1-yy).^((df2-d)/2-1).*((df2-d)/2).*yy/a;
      mu(1)=exp(gammaln(a+1)+gammaln((df2-d+2)/2)-gammaln((df2+2)/2)+a*log(df2/2));
   else
      % Density of log(exp(1)^(d/2)):
      yy=exp(y./a);   
      f(:,1)=exp(-yy).*yy/a;
      mu(1)=exp(gammaln(a+1));
   end
   
   nuv=[];
   aav=[];
   if df2<Inf
      nuv=df2+2-(1:d);
      aav=[repmat(-1/2,1,d)]; 
      for k=1:2
          if D(k)>0;
              nuv=[df1+df2-D(k) nuv];
              aav=[D(k)/2 aav];
          end;
      end;
   end
   
   for k=1:2
       if dfw1(k)<Inf & D(k)>0
           if dfw1(k)>df_limit
               nuv=[nuv dfw1(k)-dfw1(k)/dfw2(k)-(0:(D(k)-1))];
           else
               nuv=[nuv repmat(dfw1(k)-dfw1(k)/dfw2(k),1,D(k))];
           end
           aav=[aav repmat(1/2,1,D(k))];
       end
       if dfw2(k)<Inf
           nuv=[nuv dfw2(k)];
           aav=[aav -D(k)/2];
       end
   end
   
   for i=1:(numrv-1)
      nu=nuv(i);
      aa=aav(i);
      yy=y/aa+log(nu);
      % Density of log((chi^2_nu/nu)^aa):
      f(:,i+1)=exp(nu/2*yy-exp(yy)/2-(nu/2)*log(2)-gammaln(nu/2))/abs(aa);
      mu(i+1)=exp(gammaln(nu/2+aa)-gammaln(nu/2)-aa*log(nu/2));
   end
   % Check: plot(y,f); sum(f*dy,1) should be 1
      
   omega=2*pi*((1:ny)'-1)/ny/dy;
   shift=complex(cos(-b1*omega),sin(-b1*omega))*dy;
   prodfft=prod(fft(f),2).*shift.^(numrv-1);
   % Density of Y=log(B^(d/2)*U^(d/2)/sqrt(det(Q))):
   ff=real(ifft(prodfft));
   % Check: plot(y,ff); sum(ff*dy) should be 1
   mu0=prod(mu);
   % Check: plot(y,ff.*exp(y)); sum(ff.*exp(y)*dy.*(y<10)) should equal mu0   
   
   alpha=p/rhoD/mu0*fwhm(1)^D(1)*fwhm(2)^D(2)/(4*log(2))^(d/2);
   
   % Integrate the density to get the p-value for one cluster: 
   pS=cumsum(ff(ny:-1:1))*dy;
   pS=pS(ny:-1:1);
   % The number of clusters is Poisson with mean EL:
   pSmax=1-exp(-pS*EL);
   
   if p_val_extent(1) <= tlim 
      yval=minterp1(-pSmax,y,-p_val_extent);
      % Spatial extent is alpha*exp(Y) -dy/2 correction for mid-point rule:
      extent_threshold=alpha*exp(yval-dy/2);
      % For a single cluster:
      yval=minterp1(-pS,y,-p_val_extent);
      extent_threshold_1=alpha*exp(yval-dy/2);
      if length(p_val_extent)<=nprint
         extent_threshold
         extent_threshold_1
      end
   else
      % p_val_extent is now treated as a spatial extent:
      logpval=log(p_val_extent/alpha+(p_val_extent<=0))+dy/2;
      P_val_extent=interp1(y,pSmax,logpval);
      extent_threshold=P_val_extent.*(p_val_extent>0)+(p_val_extent<=0);
      % For a single cluster:
      P_val_extent_1=interp1(y,pS,logpval);
      extent_threshold_1=P_val_extent_1.*(p_val_extent>0)+(p_val_extent<=0);
      if length(p_val_extent)<=nprint
         P_val_extent
         P_val_extent_1
      end
   end
   
end
   
return

function x=gammalni(n);
i=find(n>=0);
x=Inf+n;
if ~isempty(i)
   x(i)=gammaln(n(i));
end
return

function iy=minterp1(x,y,ix);
% interpolates only the monotonically increasing values of x at ix
n=length(x);
mx=x(1);
my=y(1);
xx=x(1);
for i=2:n
   if x(i)>xx
      xx=x(i);
      mx=[mx xx];
      my=[my y(i)];
   end
end
iy=interp1(mx,my,ix);
return
