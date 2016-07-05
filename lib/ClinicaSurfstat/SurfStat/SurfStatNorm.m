function [ Y, Yav ] = SurfStatNorm( Y, mask, subdiv );

%Normalizes by subtracting the global mean, or dividing by it.
%
% Usage: [ Y, Yav ] = SurfStatNorm( Y [,mask [,subdiv ] ] );
%
% Y      = n x v matrix or n x v x k array of data, v=#vertices,
%          or memory map of same. 
% mask   = 1 x v logical vector, 1=inside, 0=outside,
%        = logical(ones(1,v)) by default.
% subdiv = 's' for Y=Y-Yav (default) or 'd' for Y=Y/Yav.
%
% Y   = n x v matrix or n x v x k array of normalized data, or memory map
%       of same. Note that if Y is a memory map, then Y is always
%       overwritten with the normalized data.
% Yav = n x 1 vector or n x k matrix of mean(Y in mask).

if isnumeric(Y)
    [n,v,k]=size(Y);
else
    Ym=Y;
    s=Ym.Format{2};
    if length(s)==2
        s=s([2 1]);
        k=1;
    else
        s=s([3 1 2]);
        k=s(3);
    end
    n=s(1);
    v=s(2);
end    
    
if nargin<2 | isempty(mask)
    mask=logical(ones(1,v));
end
if nargin<3
    subdiv = 's';
end

if isnumeric(Y)
    Yav=squeeze(mean(double(Y(:,mask,:)),2));
    for i=1:n
        for j=1:k
            if subdiv(1) == 's'
                Y(i,:,j)=Y(i,:,j)-Yav(i,j);
            else
                Y(i,:,j)=Y(i,:,j)/Yav(i,j);
            end
        end
    end
else
    Yav=zeros(n,k);
    for i=1:n
        for j=1:k
            if length(s)==2
                Y=Ym.Data(1).Data(:,i);
            else
                Y=Ym.Data(1).Data(:,j,i);
            end
            Yav(i,j)=mean(double(Y));
            if subdiv(1) == 's'
                Y=Y-Yav(i,j);
            else
                Y=Y/Yav(i,j);
            end
            if length(s)==2
                Ym.Data(1).Data(:,i)=Y;
            else
                Ym.Data(1).Data(:,j,i)=Y;
            end            
        end
    end
    Y=Ym;
end

return
end
    
    
    

