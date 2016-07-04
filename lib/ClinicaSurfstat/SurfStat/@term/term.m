function t = term( x, str );

%Makes a vector, matrix or structure into a term in a linear model.
%
% Usage: t = term( x [, str] );
%
% Internally a term consists of a cell array of strings for the names of
% the variables in the term, and a matrix with one column for each name,
% which can be accessed by char and double (see below). 
% 
% If x is a matrix of numbers, t has one variable for each column of x. 
% If x is a cell array of strings or a matrix whose rows are strings, t  
%      is a set of indicator variables for the unique strings in x. 
% If x is a structure, t has one variable for eaxh field of x. 
% If x is a scalar, t is the constant term. It is expanded to the length
%       of the other term in a binary operator.
% If x is a term, t=x. With no arguments, t is the empty term. 
% 
% str is a cell array of strings for the names of the variables. If absent,
% the names for the four cases above are either 'x' (followed by 1,2,... if
% x is a matrix), the unique strings in x, the fields of x, num2str(x) if x
% is a scalar, or '?' if x is an algebraic expression.
%
% Term t can be subscripted and returns a numeric vector or matrix, e.g.
%    t.f      = variable with name 'f'.
%    t(1,3:5) = matrix whose columns are variables 1,3,4,5.
%
% The following operators are overloaded for terms t1 and t2:
%    t1 + t2 = {variables in t1} union {variables in t2}.
%    t1 - t2 = {variables in t1} intersection complement {variables in t2}.
%    t1 * t2 = sum of the element-wise product of each variable in t1 with
%              each variable in t2, and corresponds to the interaction 
%              between t1 and t2.
%    t ^ k   = product of t with itself, k times. 
%
% Algebra: commutative, associative and distributive rules apply to terms:
%    a + b = b + a
%    a * b = b * a
%    (a + b) + c = a + (b + c)
%    (a * b) * c = a * (b * c)
%    (a + b) * c = a * c + b * c
%    a + 0 = a
%    a * 1 = a
%    a * 0 = 0
% However note that 
%    a + b - c =/= a + (b - c) =/= a - c + b
% 
% The following functions are overloaded for term t:
%    char(t)         = cell array of the names of the variables.
%    double(t)       = matrix whose columns are the variables in t, i.e.
%                      the design matrix of the linear model.
%    size(t [, dim]) = size(double(t) [, dim]).
%    isempty(t)      = 1 if t is empty and 0 if not.

if nargin == 0 
    t.names=[];
    t.matrix=[];
    t = class(t,'term');
elseif isa(x,'term')
    t = x;
elseif iscellstr(x)
    u=unique(x);
    k=length(u);
    t.names=cell(1,k);
    t.matrix=zeros(length(x),k);
    for i=1:k
        % t.names{i}=[inputname(1) '.' u{i}];
        t.names(i)=u(i);
        t.matrix(:,i)=ismember(x,u(i));
    end
    t = class(t,'term');
elseif ischar(x) && ndims(x)==2
    t.names=unique(x,'rows');
    t.matrix=zeros(length(x),length(t.names));
    for i=1:length(t.names)
        t.matrix(:,i)=ismember(x,t.names(i,:),'rows');
    end
    t = class(t,'term');
elseif isnumeric(x) && numel(x)>1
    if nargin<2
        str={inputname(1)};
%         str={'?'};
    end
    if isempty(str{1}) % changed from isempty(str) because a cell array with an empty string inside is not empty
        str={'?'};
    end
    if ischar(str)
        str={str};
    end
    if length(str)==size(x,2)
        t.names=str;
    else
        for k=1:size(x,2)
            t.names{k}=[char(str) num2str(k)];
        end
    end
    t.matrix=double(x);
    t = class(t,'term');
elseif isnumeric(x) && numel(x)==1
    t.names={num2str(x)};
    t.matrix=double(x);
    t = class(t,'term'); 
elseif isstruct(x)
    t.names=fieldnames(x)';
    t.matrix=[];
    for i=1:length(t.names)
        t.matrix=[t.matrix double(getfield(x,t.names{i}))];
    end
    t = class(t,'term'); 
end