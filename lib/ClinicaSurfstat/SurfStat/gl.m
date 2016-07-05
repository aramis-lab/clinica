function f = gl( varargin );

%Generates levels for a factor.
%
% Usage: f = gl( levels, m, n ); 
%
% levels = a single scalar, p, or
%        = a cell array of p strings {'F1', 'F2', ..., 'Fp'} for the level
%        names, or
%        = a string 'F', followed by a single scalar p, so that the call is 
%        gl( 'F', p, g, n ); the level names are 'F1', 'F2', ..., 'Fp'.
% m      = number or repetitions of each level.
% n      = length of the factor (number of observations).
%
% f =  {'F1', 'F1', ..., 'F1' (repeated m times), 
%       'F2', 'F2', ..., 'F2' (repeated m times),  
%       ...
%       'Fp', 'Fp', ..., 'Fp' (repeated m times), 
%       then all this is repeated until n observations are reached}.
% e.g. 
% gl( 3, 2, 13 ) gives
% {'F1','F1','F2','F2','F3','F3','F1','F1','F2','F2','F3','F3','F1'}
%
% gl( 's', 3, 2, 13 ) gives
% {'s1','s1','s2','s2','s3','s3','s1','s1','s2','s2','s3','s3','s1'}
%
% gl( {'a','bb','c'}, 2, 13 ) gives
% { 'a', 'a','bb','bb', 'c', 'c', 'a', 'a','bb','bb', 'c', 'c', 'a'}

if isnumeric(varargin{1})
    p=varargin{1};
    m=varargin{2};
    n=varargin{3};
    for i=1:p
        levelnames{i}=['F' num2str(i)];
    end
end
if iscell(varargin{1})
    levelnames=varargin{1};
    p=length(levelnames);
    m=varargin{2};
    n=varargin{3};
end
if isstr(varargin{1})
    p=varargin{2};
    m=varargin{3};
    n=varargin{4};
    for i=1:p
        levelnames{i}=[varargin{1} num2str(i)];
    end
end

[a,b,c]=ind2sub([m p ceil(n/(m*p))],1:n);
f=levelnames(b)';

return
end
    


