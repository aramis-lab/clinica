function s=mpower(t,p)
t=term(t,inputname(1));
if p>=2
    s=mtimes(t,mpower(t,p-1));
else 
    s=t;
end
    