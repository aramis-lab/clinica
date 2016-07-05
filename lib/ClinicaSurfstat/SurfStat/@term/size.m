function s=size(t,dim)
if nargin==1
    s=size(t.matrix);
else
    s=size(t.matrix,dim);
end