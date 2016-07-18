function s=size(m,dim)
s=[size(m.mean) size(m.variance)];
if nargin==2
    s=s(dim);
end