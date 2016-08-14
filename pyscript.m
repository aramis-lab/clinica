fprintf(1,'Executing %s at %s:\n',mfilename,datestr(now));
ver,
try,
addpath('/home/guest/.local/lib/python2.7/site-packages/Clinica-0.1.0-py2.7.egg/clinica/lib/clinicasurfstat');

    clinicasurfstat('/home/guest/HAO/Clinica/clinica/examples/external-data/clinica_surfstat', '/tmp/tmpSaZDkm', '1 + Label + Gender + Age', 'Label', '/home/guest/HAO/Clinica/clinica/examples/external-data/clinica_surfstat/csv_file/template.csv', '%s %s %s %f');
    
,catch ME,
fprintf(2,'MATLAB code threw an exception:\n');
fprintf(2,'%s\n',ME.message);
if length(ME.stack) ~= 0, fprintf(2,'File:%s\nName:%s\nLine:%d\n',ME.stack.file,ME.stack.name,ME.stack.line);, end;
end;