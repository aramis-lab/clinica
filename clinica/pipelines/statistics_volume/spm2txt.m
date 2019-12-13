function spm2txt(filename, xSPM)
%o=spm2txt(filename)
%
%Export the current SPM table to a .txt file FILENAME. Import to your
%spreadsheet engine as a tab separated data, and voila.
%
%Selim Onat (2018)

%%
t           = spm_list('table',xSPM);
for C = 1:size(t.dat,2)
    o(:,C)  = cellfun(@fun,t.dat(:,C), 'UniformOutput',false);%fun defined below
end
t.hdr(1,[2 4 5 6 7 8 10 11]) = {'' '' '' '' '' '' '' ''};%remove unnecessary headers
o  = o';
%% open, write, close file
f  = fopen(filename,'w');
fprintf(f,['Statistics: %-11s\t\n'],t.tit);
fprintf(f,[repmat('%-11s\t',1,12) '\n'],t.hdr{1,:});
fprintf(f,[repmat('%-11s\t',1,12) '\n'],t.hdr{2,:});
fprintf(f,[repmat('=',1,12*12) '\n']);
fprintf(f,[repmat('%-11s\t',1,12) '\n'],o{:});
fprintf(f,'\n\n\n%s\n',t.str);
fprintf(f,[repmat('=',1,12*12) '\n']);
for foot = 1:size(t.ftr,1)
    fprintf(f,[t.ftr{foot,1} '\n'],t.ftr{foot,2});
end   
fclose(f);
%%
    function o=fun(x)
        %transform empty cells to ''
        if isempty(x)
            o='';
        else
            o=sprintf(t.fmt{C},x);
        end
    end

end