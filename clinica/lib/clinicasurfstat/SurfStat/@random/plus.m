function s=plus(m1,m2)

if (~isa(m1,'term') && ~isa(m1,'random') && numel(m1)>1) || ...
   (~isa(m2,'term') && ~isa(m2,'random') && numel(m2)>1)
    warning('If you don''t convert vectors to terms you can get unexpected results :-(') 
end
if ~isa(m1,'random')
    m1=random([],m1,[],inputname(1));
end
if ~isa(m2,'random')
    m2=random([],m2,[],inputname(2));
end

if size(m1,3)==1
    v=eye(max(size(m2,1),sqrt(size(m2,3))));   
    m1.variance=term(v(:),'I');
end
if size(m2,3)==1
    v=eye(max(size(m1,1),sqrt(size(m1,3))));
    m2.variance=term(v(:),'I');
end

s.mean=m1.mean+m2.mean;
s.variance=m1.variance+m2.variance;
s=class(s,'random');