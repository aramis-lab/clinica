function [surf,ab]=SurfStatReadSurf(filenames,ab,numfields,dirname,maxmem);

%Reads coordinates and triangles from an array of .obj or FreeSurfer files. 
%
% Usage: [ surf, ab ] = SurfStatReadSurf( filenames [,ab [,numfields ...
%                                [, dirname [, maxmem ] ] ] ] );
% 
% filenames = .obj or FS file name (n=1) or n x k cell array of file names.
% ab        = 'a' for ASCII or 'b' for binary. If it doesn't work it 
%             will try the other. Default is 'a'.
% numfields = number of fields to read, in the order below, default 2.
% dirname   = name of a directory where you have write permission. This
%             is only needed if the data is memory mapped. The data is
%             memory mapped only if it exceeds maxmem Mb as 4 byte reals.
%             SurfStatReadSurf then writes the data as 4 byte reals to a
%             temporary file in dirname, then memory maps this file to the
%             output data. If dirname does not exist, SurfStatReadVol will
%             create it. The default is
%           = (filenames{1,1} directory)/SurfStat, so you can ignore this
%             parameter if you have write permission in the filenames{1,1} 
%             directory.
% maxmem    = memory limit in Mb. The default is 64. 
%
% surf.tri   = t x 3 matrix of triangle indices, 1-based, t=#triangles.
% surf.coord = 3 x v matrix of coordinates, v=#vertices, if n=1, or
%              n x v x 3 array if n>1, or memory map of same. Data from the
%              k files are concatenated. Note that the mapped file is not
%              deleted after you quit MATLAB.
% ab         = whichever was successful.
% The coordinates and triangle indices of the k files are concatenated. 
% If n=1, surf.coord is double precision; 
% if n>1, surf.coord is single precision.
% surf.tri is int32.

global MAXMEM
if nargin<5
    if isempty(MAXMEM)
        maxmem=64;
    else
        maxmem=MAXMEM;
    end
end
maxmem=maxmem*2^20;

if nargin<2 | isempty(ab)
    ab='a';
end
if nargin<3 | isempty(numfields)
    numfields=2;
end
numfields1=(numfields-1)*3+1;

if isstr(filenames)
    sf=filenames;
    filenames=cell(1,1);
    filenames(1)={sf};
end

[n,k]=size(filenames);
if n>1
    fprintf(1,'%s',[num2str(n) ' x ' num2str(k) ' files to read, % remaining: 100 ']);
end
c=[];
vs=zeros(1,k);
if numfields==2
    surf.tri=[];
end
for j=1:k
    [s,ab]=SurfStatReadSurf1(filenames{1,j},ab,numfields1);
    if numfields==2
        surf.tri=[surf.tri; int32(s.tri)+size(c,2)];
    end
    c=[c s.coord];
    vs(j)=size(s.coord,2);
end

if n==1
    surf.coord=c;
else
    n10=floor(n/10);
    v=size(c,2);
    isnum=(n*v*3*4<=maxmem);
    if isnum
        surf.coord=zeros(n,v,3,'single');
        surf.coord(1,:,:)=single(c)';
    else
        if nargin<4 | isempty(dirname)
            [PATHSTR,NAME,EXT]=fileparts(filenames{1,1});
            dirname=fullfile(PATHSTR,'SurfStat');
        end
        if ~exist(dirname,'dir')
            [SUCCESS,MESSAGE,MESSAGEID]=mkdir(dirname);
            if ~SUCCESS
                error(MESSAGEID,['Tried to make directory ' dirname ' for memory mapping. \n',...
                    'Please specify the name of a directory where you have write permission \n',...
                    'as the fourth parameter of SurfStatReadSurf - see the help :-) Error: \n',...
                    MESSAGE]);
            end
        end
        [PATHSTR,NAME,EXT]=fileparts(tempname);
        Filename=fullfile(dirname,NAME);
        fid=fopen(Filename,'wb');
        fwrite(fid,c','single');
    end
    v2=cumsum(vs);
    v1=v2-vs+1;
    for i=2:n
        if rem(i,n10)==0 
            fprintf(1,'%s',[num2str(round(100*(1-i/n))) ' ']);
        end
        for j=1:k
            [s,ab]=SurfStatReadSurf1(filenames{i,j},ab,1);
            if isnum
                surf.coord(i,v1(j):v2(j),:)=single(s.coord)';
            else
                c(:,v1(j):v2(j))=single(s.coord);
            end
        end
        if ~isnum
            fwrite(fid,c','single');
        end
    end
    fprintf(1,'%s\n','Done');
    if ~isnum
        fclose(fid);
        surf.coord=memmapfile(Filename,'Format',{'single' [v 3 n] 'Data'},'Writable',true);
        surf.coord
    end 
end

return
end

