function s=mpower(m,p)
if p>=2
    s=mtimes(m,mpower(m,p-1));
else 
    s=m;
end