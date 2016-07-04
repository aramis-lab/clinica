function filenames = SurfStatListDir( d, exclude );

%Lists all the file names in a directory.
%
% Usage: filenames = SurfStatListDir( dir [, exclude] );
%
% d         = directory, including wildcards such as '*'.
% exclude   = matrix with rows as strings; any file name containing
%             any of these strings will be excluded. Empty by default.
%
% filenames = structure array of file names.
%
% e.g. SurfStatListDir('c:/somedir/someplace/subject*_thickness.txt');

if nargin<2
    exclude=[];
end
path=fileparts(d);
f=dir(d);
j=0;
for i=1:length(f);
    if ~f(i).isdir 
        y=1;
        for k=1:size(exclude,1)
            y=y&isempty(strfind(f(i).name,exclude(k,:)));
        end
        if y
        j=j+1;
        filenames{j,1}=fullfile(path,f(i).name);
        end
    end
end
fprintf(1,'%s\n',['List ', num2str(j) ' file names, first one is:']);
fprintf(1,'%s\n',filenames{1});

return
end
