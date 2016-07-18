function s=mtimes(m1,m2)

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

s.mean=m1.mean*m2.mean;
s.variance=m1.variance*m2.variance;

N=char(m1.mean);
if ~isempty(N)
    X=double(m1.mean);
    X=X./repmat(max(abs(X)),size(X,1),1);
    k=length(N);
    t=term;
    for i=1:k;
        for j=1:i;
            if i==j
                v=X(:,i)*X(:,i)';
                t=t+term(v(:),N{i});
            else
                v=(X(:,i)+X(:,j))*(X(:,i)+X(:,j))'/4;
                t=t+term(v(:),['(' N{j} '+' N{i} ')']);
                v=(X(:,i)-X(:,j))*(X(:,i)-X(:,j))'/4;
                t=t+term(v(:),['(' N{j} '-' N{i} ')']);
            end
        end
    end
    s.variance=s.variance+t*m2.variance;
end

N=char(m2.mean);
if ~isempty(N)
    X=double(m2.mean);
    X=X./repmat(max(abs(X)),size(X,1),1);
    k=length(N);
    t=term;
    for i=1:k;
        for j=1:i;
            if i==j
                v=X(:,i)*X(:,i)';
                t=t+term(v(:),N{i});
            else
                v=(X(:,i)+X(:,j))*(X(:,i)+X(:,j))'/4;
                t=t+term(v(:),['(' N{j} '+' N{i} ')']);
                v=(X(:,i)-X(:,j))*(X(:,i)-X(:,j))'/4;
                t=t+term(v(:),['(' N{j} '-' N{i} ')']);
            end
        end
    end
    s.variance=s.variance+m1.variance*t;
end

s=class(s,'random');

return
end