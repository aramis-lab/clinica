function display(model)
display(model.mean);
names=char(model.variance);
v=double(model.variance);
[n2,l]=size(v);
n=sqrt(n2);
for k=1:l
    disp(names{k});
%     d=[repmat('| ',n,1) num2str(reshape(v(:,k),n,n)) repmat(' |',n,1)];
%     disp(['+' repmat('-',1,size(d,2)-2) '+']);
%     disp(d);
%     disp(['+' repmat('-',1,size(d,2)-2) '+']);
    d=num2str(reshape(v(:,k),n,n));
    disp(repmat('v',1,size(d,2))); 
    disp(d);
    disp(' ');
end
