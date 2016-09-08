function [out,LINE] = go_to_key_word(fid, keyWord)
  
% Usage: [out, line] = go_to_key_word( fid, keyWord)
%
% Input: 
%  - fid: file descriptor, as obtained with fopen();
%  - keyWord: requested key word;
% 
% Output:
%  - out:  key word actually found (-1 if not found);
%  - line: number of the line in the file where the key word was found;
%
% Copyright Pierre Fillard - INRIA Sophia Antipolis, FRANCE - 2005.
%
  
frewind(fid);
LINE = 1;
newL = lower( fgetl(fid) );

while( ~strncmp(newL, lower(keyWord), length(keyWord)) )
    newL = lower( fgetl(fid) );
    if( feof(fid) )
      out = -1;
      return;
    end
    LINE = LINE+1;
end

out = newL;