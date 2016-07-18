function r=subsref(t,s)
switch s.type
    case '()'
        cols=s.subs{:};
    case '.'
        cols=find(ismember(t.names,s.subs));
end
r=t.matrix(:,cols);        


