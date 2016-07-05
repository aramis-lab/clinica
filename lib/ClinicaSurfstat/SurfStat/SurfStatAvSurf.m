function surf = SurfStatAvSurf( filenames, fun );

%Average, minimum or maximum of surfaces.
%
% Usage: surf = SurfStatAvSurf( filenames [, fun] );
%
% filenames = n x k cell array of names of files.
% fun       = function handle to apply to two surfaces, e.g. 
%           = @plus (default) will give the average of the surfaces,
%           = @min or @max will give the min or max, respectively. 
%
% surf.coord = 3 x v matrix of average coordinates, v=#vertices.
% surf.tri   = t x 3 matrix of triangle indices, 1-based, t=#triangles.
% The coordinates and triangle indices of the k files are concatenated. 

if nargin<2
    fun=@plus;
end

[n,k]=size(filenames);
fprintf(1,'%s',[num2str(n) ' x ' num2str(k) ' files to read, % remaining: 100 ']);
n10=floor(n/10);
ab='a';
for i=1:n
    if rem(i,n10)==0
        fprintf(1,'%s',[num2str(round(100*(1-i/n))) ' ']);
    end
    if i==1
        [s,ab]=SurfStatReadSurf(filenames(i,:),ab,2);
        surf.tri=s.tri;
        surf.coord=double(s.coord);  
        m=1;
% passing ab speeds up reading if the next file has the same format. 
    else
        [s,ab]=SurfStatReadSurf(filenames(i,:),ab,1);
        surf.coord=fun(surf.coord,double(s.coord));
        m=fun(m,1);
    end
end
surf.coord=surf.coord/m;
fprintf(1,'%s\n','Done');  

return
end
