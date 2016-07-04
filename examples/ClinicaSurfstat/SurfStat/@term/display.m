function display(t)
disp(' ');
if isempty(t.matrix)
    disp('[]');
else
    d=[];
    n=size(t.matrix,1);
    for k=1:length(t.names)
        d=[d repmat('  ',n+1,1) strvcat(t.names{k},num2str(t.matrix(:,k)))];
    end
    d=[d(1,:); repmat('-',1,size(d,2)); d((1:n)+1,:)];
    disp(d);
end
disp(' ');