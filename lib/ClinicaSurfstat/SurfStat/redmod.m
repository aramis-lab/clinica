function mr = redmod( m, meanorvariance );

%Reduces a linear mixed effects model by removing redundant variables. 
%
% Usage: mr = redmod( m [,meanorvariance] );
%
% m = model, either term or random or anything that can be so converted.
% meanorvariance = 'm' to reduce just the mean, 'v' to reduce just the
%                  variance, or both if absent. 
%
% mr = reduced model. Variables with a non-zero coefficient in the null  
%      space of the design matrix, and with the largest number of non-zero 
%      entries are eliminated first, until the design matrix has full rank. 

if ~isa(m,'term') && ~isa(m,'random') && numel(m)>1
    warning('If you don''t convert vectors to terms you can get unexpected results :-(') 
end
if ~isa(m,'random')
    m=random([],m,[],inputname(1));
end
[X,V]=double(m);
[C,D]=char(m);

fprintf(1,'%s','#variables remaining in model for:');
if nargin<2 | lower(meanorvariance(1))=='m'
    p=size(X,2);
    fprintf(1,'%s',[' Mean ' num2str(p)]);
    nX=null(X,'r')';
    while ~isempty(nX)
        nn=size(nX,1);
        a=repmat(mean(abs(X))./max(abs(X)),nn,1).*(abs(nX)>eps);
        [y,i]=max(a');
        keep=ones(1,p)>0;
        keep(i)=0;
        X=X(:,keep);
        C=C(keep);
        p=size(X,2);
        fprintf(1,'%s',[' ' num2str(p)]);
        nX=null(X,'r')';
    end
    fprintf(1,'%s',' Done.');
end

if nargin<2 | lower(meanorvariance(1))=='v'
    q=size(V,2);
    fprintf(1,'%s',[' Variance: ' num2str(q)]);
    nV=null(V,'r')';
    while ~isempty(nV)
        nn=size(nV,1);
        a=repmat(mean(abs(V))./max(abs(V)),nn,1).*(abs(nV)>eps);
        [y,i]=max(a');
        keep=ones(1,q)>0;
        keep(i)=0;
        V=V(:,keep);
        D=D(keep);
        q=size(V,2);
        fprintf(1,'%s',[' ' num2str(q)]);
        nV=null(V,'r')';
    end
    fprintf(1,'%s',' Done.');
end

if isempty(V)
    mr=term(X,C);
else
    mr=random(term(V,D),term(X,C),[],[],1);
end
fprintf(1,'\n');

return
end
