function data = SurfStatReadData( filenames, dirname, maxmem )

%Reads data (e.g. thickness) from an array of .txt or FreeSurfer files. 
%
% Usage: data = SurfStatReadData( filenames [, dirname [, maxmem ]] );
%
% filenames = .txt, .mgh or FS file name (n=1) or n x k cell array of file 
%             names. If the .mgh file has data from n surfaces, then
%             filenames must be a file name or 1 x k cell array of file 
%             names.
% dirname   = name of a directory where you have write permission. This
%             is only needed if the data is memory mapped. The data is
%             memory mapped only if it exceeds maxmem Mb as 4 byte reals.
%             SurfStatReadData then writes the data as 4 byte reals to a
%             temporary file in dirname, then memory maps this file to the
%             output data. If dirname does not exist, SurfStatReadVol will
%             create it. The default is
%           = (filenames{1,1} directory)/SurfStat, so you can ignore this
%             parameter if you have write permission in the filenames{1,1} 
%             directory.
% maxmem    = memory limit in Mb. The default is the gloabl MAXMEM, or 64. 
%
% data = n x v matrix of data, v=#vertices, or memory map of same. Data
%        from the k files are concatenated. Note that the mapped file is
%        not deleted after you quit MATLAB.
% If n=1, data is double precision; 
% if n>1, data is single precision.

global MAXMEM
% If the global variable doesn't exist the first time you issue
%the global statement, it will be initialized to the empty matrix.
if nargin<3
    if isempty(MAXMEM)
        maxmem=64;
    else
        maxmem=MAXMEM;
    end
end
maxmem=maxmem*2^20;

if isstr(filenames)
    sf=filenames;
    filenames=cell(1,1);
    filenames(1)={sf};
end

[n,k]=size(filenames);
if n>1
    fprintf(1,'%s',[num2str(n) ' x ' num2str(k) ' files to read, % remaining: 100 ']);
end
data=[];
vs=zeros(1,k);
for j=1:k
    data1=SurfStatReadData1(filenames{1,j});% this is to read the data from .txt/.thickness/.mgh 
    if size(data1,1)==1
        data=[data data1];
    else
        data=[data single(data1)];
    end        
    vs(j)=size(data1,2);
end

if n>1
    n10=floor(n/10);
    data1=data;
    v=size(data1,2);
    isnum=(n*v*4<=maxmem);
    if isnum
        data=zeros(n,v,'single');
        data(1,:)=single(data1);
    else
        if nargin<2 | isempty(dirname)
            [PATHSTR,NAME,EXT]=fileparts(filenames{1,1});
            dirname=fullfile(PATHSTR,'SurfStat');
        end
        if ~exist(dirname,'dir')
            [SUCCESS,MESSAGE,MESSAGEID]=mkdir(dirname);
            if ~SUCCESS
                error(MESSAGEID,['Tried to make directory ' dirname ' for memory mapping. \n',...
                    'Please specify the name of a directory where you have write permission \n',...
                    'as the second parameter of SurfStatReadData - see the help :-) Error: \n',...
                    MESSAGE]);
            end
        end
        [PATHSTR,NAME,EXT]=fileparts(tempname);
        Filename=fullfile(dirname,NAME);
        fid=fopen(Filename,'wb');
        fwrite(fid,data1,'single');
    end
    v2=cumsum(vs);
    v1=v2-vs+1;
    for i=2:n
        if rem(i,n10)==0 
            fprintf(1,'%s',[num2str(round(100*(1-i/n))) ' ']);
        end
        for j=1:k
            data1=SurfStatReadData1(filenames{i,j});
            if isnum
                data(i,v1(j):v2(j))=single(data1);
            else
                fwrite(fid,data1,'single');
            end
        end
    end
    fprintf(1,'%s\n','Done');
    if ~isnum
        fclose(fid);
        data=memmapfile(Filename,'Format',{'single' [v n] 'Data'},'Writable',true)
    end
end
    
return
end


