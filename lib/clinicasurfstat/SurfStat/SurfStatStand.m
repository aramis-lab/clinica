function [ Y, Ym ] = SurfStatStand( Y, mask, subtractordivide);

%Standardizes by subtracting the global mean, or dividing by it.
%
% Usage: [ Y, Ym ] = SurfStatStand( Y [,mask [,subtractordivide] ] );
%
% Y                = n x v matrix of data, v=#vertices.
% mask             = 1 x v, 1=inside, 0=outside, v=#vertices.
%                  = ones(1,v) by default.
% subtractordivide = 's' for Y=Y-Ym (default) or 'd' for Y=(Y/Ym-1)*100.
%
% Y  = n x v matrix of standardized data.
% Ym = n x 1 vector of mean(Y in mask).

if nargin<2 | isempty(mask)
    mask=logical(ones(1,size(Y,2)));
end
if nargin<3
    subtractordivide = 's';
end
Ym=mean(Y(:,mask),2);
for i=1:size(Y,1)
    if subtractordivide(1) == 's'
        Y(i,:)=Y(i,:)-Ym(i);
    else
        Y(i,:)=(Y(i,:)/Ym(i)-1)*100;
    end
end

return
end
    
    
    

