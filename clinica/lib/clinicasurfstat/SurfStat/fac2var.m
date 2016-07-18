function x = fac2var( f, v )

%Converts a cell array of strings or term of this to a numeric variable. 
%
% Usage: x = fac2var( f [,v] );
%
% f = n x 1 cell array of strings or term of this.
% v = 1 x k numeric vector of values for the k unique strings in f, in 
%     ascending order, =1:k if absent. 
%
% x = n x 1 numeric vector.

if isa(f,'term')
    f=double(f)*(1:size(f,2))';
end
[u,i,j]=unique(f);
if nargin<2
    v=1:length(u);
end
x=v(j)';
