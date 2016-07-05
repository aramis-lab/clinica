function f = var2fac( x, str );

%Converts a numeric variable or term of this to a cell array of strings.
%
% Usage: f = var2fac( x [, str] );
%
% x   = n x 1 numeric vector or term of this.
% str = either a single string, or a 1 x k cell array of strings of names 
%       for the k unique values of x in ascending order. If str is absent,
%       then str='x'.
%
% f = n x 1 cell array of strings. If str is a single string, then the 
%     strings are 'str1', 'str2', .... 'strk'.
% 
% Examples: 
% gender=var2fac([repmat(1,1,64) repmat(2,1,44)],{'male'; 'female'})
% subj=var2fac(repmat(1:27,4,1),'subj')

x=double(x);
n=numel(x);
if nargin<2
    str=inputname(1);
end
if ischar(str)
    for i=1:n
        f{i}=[str num2str(x(i))];
    end
    f=f';
else
    [u,i,j]=unique(x);
    f=str(j);
end
    