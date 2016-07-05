function s=mtimes(t1,t2)
if (~isa(t1,'term') && numel(t1)>1) | (~isa(t2,'term') && numel(t2)>1)
    warning('If you don''t convert vectors to terms you can get unexpected results :-(') 
end
t1=term(t1,inputname(1));
t2=term(t2,inputname(2));
if isempty(t1) | isempty(t2)
    s=term;
    return
end
n1=size(t1.matrix,1);
n2=size(t2.matrix,1);
n=max(n1,n2);
if n1==1
    t1.matrix=repmat(t1.matrix,n,1);
end
if n2==1
    t2.matrix=repmat(t2.matrix,n,1);
end
k1=size(t1.matrix,2);
k2=size(t2.matrix,2);
nms=cell(1,k1*k2);
x=zeros(n,k1*k2);
for i2=1:k2
    for i1=1:k1
        i1i2=i1+(i2-1)*k1;
        if ~any(t1.matrix(:,i1)-1)
            nms(i1i2)=t2.names(i2);
            x(:,i1i2)=t2.matrix(:,i2);
        elseif ~any(t2.matrix(:,i2)-1)
            nms(i1i2)=t1.names(i1);
            x(:,i1i2)=t1.matrix(:,i1);
        else
            nms(i1i2)={[t1.names{i1} '*' t2.names{i2}]};
            x(:,i1i2)=t1.matrix(:,i1).*t2.matrix(:,i2);
        end
    end
end
m=x./repmat(sum(abs(x)),n,1);
[u,j]=unique(flipud(m'),'rows');
j=sort(k1*k2-j+1);
jj=j(any(x(:,j)));
s.names=nms(jj);
s.matrix=x(:,jj);
s=class(s,'term');
