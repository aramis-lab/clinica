function SurfStatDelete( varargin );

%Deletes variables including memory mapped files.
%
% Usage: SurfStatDelete( variable_list [, qualifiers] );
%
% variable_list = comma-delimited list of quoted strings: 'var1', 'var2',
% ..., 'varN'. You can use the wildcard character * to delete variables
% that match a pattern. For example, SurfStatDelete('A*') deletes all
% variables in the current workspace that start with A.
%
% If a variable is a memory map, then the mapped file is also deleted. 
%
% To clean up everything, use SurfStatDelete with no arguments (you will be
% prompted to confirm). It is a good idea to do this before exiting MATLAB,
% otherwise memory mapped files will not be deleted. 
%
% qualifiers = those variables in variable_list that meet all
% qualifications specified in qualifiers. See "help who" for a list.

if nargin==0
    reply=input('Are you sure you want to delete everything? Y/N [Y]: ', 's');
    if isempty(reply)
        reply='Y';
    end
    if ~strcmp(reply,'Y')
        return
    end
end
if nargin==0
    w=evalin('caller','whos');
else
    s=['''' varargin{1} ''''];
    for i=2:nargin
        s=[s ' , ''' varargin{i} ''''];
    end
    w=evalin('caller',['whos(' s ')']);
end
n=length(w);
for i=1:n
    isc=false;
    if strcmp(w(i).class,'struct')
        if evalin('caller',['isfield(' w(i).name ',''coord'')'])
            isc=evalin('caller',['isa(' w(i).name '.coord,''memmapfile'')']);
        end
    end
    if strcmp(w(i).class,'memmapfile') | isc
        if isc 
            m=evalin('caller',[w(i).name '.coord']);
        else          
            m=evalin('caller',w(i).name);
        end
        filename=m.Filename;
        evalin('caller',['clear ' w(i).name]);        
        [pathstr,name,ext]=fileparts(filename);
        filenameResid=fullfile(pathstr,['Resid' ext]);
        warning off all;
        delete(filename);
        delete(filenameResid);       
        warning on all;
        if length(dir(pathstr))==2
            rmdir(pathstr);
        end
    else
        evalin('caller',['clear ' w(i).name]);
    end
end

return
end