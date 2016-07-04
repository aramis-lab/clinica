function d = SurfStatReadVol1( file, Z, T );

%Reads a single volumetric file in MINC, ANALYZE, NIFTI or AFNI format. 
% 
% Usage: d = SurfStatReadVol1( file [, Z, T] ). 
%
% file = file name with extension .mnc, .img, .nii or .brik as above.
% Z    = vector of slices.
% T    = vector of times. 
% If Z and T are both 0, then just the header is read into d. 
% If Z and T are omitted, reads the entire volume.

if ~isstr(file)
    d = file;
    if isfield(d,'data') & nargin==3
        if length(size(d.data))==4 & Z~=0 & T~=0
            d.data=file.data(:,:,Z,T);
        end
        if length(size(d.data))==3 & length(Z)>1 & length(T)==1
            d.data=file.data(:,:,Z);
        end
        if length(size(d.data))==3 & length(Z)==1 & length(T)>1
            d.data=file.data(:,:,T);
        end
        if length(size(d.data))==3 & length(Z)==1 & length(T)==1 & Z>0 & T>0
            d.data=file.data(:,:,Z);
        end
    end
    return
end

file = deblank(file);
   
[path,name,ext]=fileparts(file);
if strcmp(ext,'.gz')
   if isunix
      unix(['gunzip ' file]);
   else
      ['Can''t gunzip on non-unix system']
      return
   end
   [path,name,ext]=fileparts(name)
end 

switch lower(ext)
case '.mnc'
   if nargin == 3
      d = fmris_read_minc(file,Z,T);
   else
      d = fmris_read_minc(file);
   end
case '.img'
   if nargin == 3
      d = fmris_read_analyze(file,Z,T);
   else
      d = fmris_read_analyze(file);
   end
   if isfield(d,'data') d.data(isnan(d.data))=0; end;
   d.origin=-d.origin.*d.vox(1:3);
case '.brik'
   if nargin == 3
      d = fmris_read_afni(file,Z,T);
   else
      d = fmris_read_afni(file);
   end
case '.nii'
   if nargin == 3
      d = fmris_read_nifti(file,Z,T);
   else
      d = fmris_read_nifti(file);
   end
otherwise
   ['Unknown file extension']
end

d.parent_file=file;

return
end



%%

function [d]=fmris_read_afni(file,Z,T)

[err, Info]=BrikInfo(file);

d.dim 	= [Info.DATASET_DIMENSIONS(1:3) Info.DATASET_RANK(2)];

if nargin == 1 
   Z = 1:d.dim(3);
   T = 1:d.dim(4);
end

if (T(1)~=0)&(Z(1)~=0)
   Opt.Slices=Z;
   Opt.Frames=T;
   [err, d.data, Info, ErrMessage] = BrikLoad(file, Opt);
   d.calib		= [min(min(min(min(d.data)))) max(max(max(max(d.data))))];
end

d.vox 					= Info.DELTA;
d.vox_units				= '';
d.vox_offset			= 0;
d.precision				= '';
d.calib_units			= '';
d.origin 				= Info.ORIGIN;
d.descrip				= '';

if isfield(Info,'WORSLEY_DF')
   df=Info.WORSLEY_DF;
else
   df=[];
end
if ~isempty(df)
   d.df=df;
end

if isfield(Info,'WORSLEY_NCONJ')
   nconj=Info.WORSLEY_NCONJ;
else
   nconj=[];
end
if ~isempty(nconj)
   d.nconj=nconj;
end

if isfield(Info,'WORSLEY_FWHM')
   fwhm=Info.WORSLEY_FWHM;
else
   fwhm=[];
end
if ~isempty(fwhm)
   d.fwhm=fwhm;
end

return;
end



%%

function	d=fmris_read_analyze(file,Z,T);

% Read in 1 4D or >1 3D Analyze format Image Volumes
%
% Usage
%
% d=Read_Analyze(file);
% Reads in all image data and attributes into the structure d.
%
% d=Read_Analyze(file,Z,T);
% Reads in chosen planes and frames and attributes into the structure d.
% Z and T are vectors e.g. Z=1:31;T=1:12; or Z=[24 26];T=[1 3 5];
%
% (c) Roger Gunn & John Aston  

% Get machine format on which the file was written:

machineformat=getmachineformat(file);

if ~isempty(machineformat)
   
   % Read Header Information
   d=Read_Analyze_Hdr(file,machineformat);
   d.file_name=file;
   
   % try to get file precision if it is unknown:
   if strcmp(d.precision,'Unknown')
      d.precision=getprecision(deblank(file),machineformat,d.dim,d.global);
   end
   
   % Read Image Data
   
   if ~isempty(d) & ~strcmp(d.precision,'Unknown') 
      
      if nargin==1 | strcmp(d.precision,'uint1') 
         
         % Read in Whole Analyze Volume
         fid = fopen(d.file_name,'r',machineformat);
         if fid > -1
            d.data=d.scale*reshape(fread(fid,prod(d.dim),d.precision),d.dim(1),d.dim(2),d.dim(3),d.dim(4));
            fclose(fid);
         else
            errordlg('Check Image File: Existence, Permissions ?','Read Error'); 
         end;
         
         if nargin==3
            if all(Z>0)&all(Z<=d.dim(3))&all(T>0)&all(T<=d.dim(4))
               d.data=d.data(:,:,Z,T);
               d.Z=Z;
               d.T=T;
            else
               errordlg('Incompatible Matrix Identifiers !','Read Error');  
            end
         end
         
      elseif nargin==3
         % Read in Chosen Planes and Frames
         if (T(1)~=0)|(Z(1)~=0)
            
            if all(Z>0)&all(Z<=d.dim(3))&all(T>0)&all(T<=d.dim(4))
               
               fid = fopen(d.file_name,'r',machineformat);
               if fid > -1
                  d.data=zeros(d.dim(1),d.dim(2),length(Z),length(T));
                  for t=1:length(T)
                     for z=1:length(Z)
                        status=fseek(fid,d.hdr.byte*((T(t)-1)*prod(d.dim(1:3))+(Z(z)-1)*prod(d.dim(1:2))),'bof');
                        d.data(:,:,z,t)=d.scale*fread(fid,[d.dim(1) d.dim(2)],d.precision);
                     end
                  end
                  d.Z=Z;
                  d.T=T;
                  fclose(fid);
               else
                  errordlg('Check Image File: Existence, Permissions ?','Read Error'); 
               end;
               
            else
               errordlg('Incompatible Matrix Identifiers !','Read Error'); 
            end;
            
         end;
      else
         errordlg('Unusual Number of Arguments','Read Error');
      end;
   else
      if strcmp(d.precision,'Unknown');
         errordlg('Unknown Data Type (Precision?)','Read Error');
      end
   end;
else
   errordlg('Unknown Machine Format','Read Error'); 
end

% if there is no slice thickness, set it to 6mm:

if d.vox(3)==0
   d.vox(3)=6;
end

return;
end

%%

function machineformat=getmachineformat(file);

% Get machine format by reading the d.hdr.dim(1) attribute and 
% making sure it is 1, 2, 3 or 4. 

machineformat=[];
for mf='nlbdgcas'
   fid = fopen([file(1:(length(file)-3)) 'hdr'],'r',mf);
   if fid > -1
      fseek(fid,40,'bof');
      if any(fread(fid,1,'int16')==1:4)
         machineformat=mf;
         fclose(fid);
         break
      else
         fclose(fid);
      end
   else
      errordlg('Check Header File: Existence, Permissions ?','Read Error'); 
      break
   end
end

return
end

%%

function precision=getprecision(file,machineformat,dim,range);

% Get precision by reading a value from the middle of the .img file and 
% making sure it is within the global attribute

precisions=['int8   '   
        'int16  '  
        'int32  '  
        'int64  '   
        'uint8  '   
        'uint16 '  
        'uint32 '  
        'uint64 ' 
        'single ' 
        'float32' 
        'double ' 
        'float64' ];
nbytes=[1 2 4 8 1 2 4 8 4 4 8 8];
middle_vol=dim(1)*dim(2)*floor(dim(3)/2)+dim(1)*round(dim(2)/2)+round(dim(1)/2);
h=dir(file);
n=ceil(h.bytes/prod(dim));

fid = fopen(file,'r',machineformat);
if fid > -1
   for i=1:size(precisions,1)
      if nbytes(i)==n
         status=fseek(fid,middle_vol*n,'bof');
         if status==0
            precision=deblank(precisions(i,:));
            x=fread(fid,10,precision);
            if all(range(1)<=x) & all(x<=range(2))
               return
            end
         end
      end
   end
end
errordlg('Check Header File: Existence, Permissions ?','Read Error'); 
precision='Unknown';

return
end

%%

function d=Read_Analyze_Hdr(file,machineformat);

% Read Analyze Header information into the structure d
% Adapted from John Ashburners spm_hread.m

fid  = fopen([file(1:(length(file)-3)) 'hdr'],'r',machineformat);

if fid > -1
   
   % read (struct) header_key
   %---------------------------------------------------------------------------
   fseek(fid,0,'bof');
   
   d.hdr.sizeof_hdr 		= fread(fid,1,'int32');
   d.hdr.data_type  		= deblank(setstr(fread(fid,10,'char'))');
   d.hdr.db_name    		= deblank(setstr(fread(fid,18,'char'))');
   d.hdr.extents    		= fread(fid,1,'int32');
   d.hdr.session_error   	= fread(fid,1,'int16');
   d.hdr.regular    		= deblank(setstr(fread(fid,1,'char'))');
   d.hdr.hkey_un0    		= deblank(setstr(fread(fid,1,'char'))');
   
   
   
   % read (struct) image_dimension
   %---------------------------------------------------------------------------
   fseek(fid,40,'bof');
   
   d.hdr.dim    			= fread(fid,8,'int16');

   d.hdr.vox_units    		= deblank(setstr(fread(fid,4,'char'))');
   d.hdr.cal_units    		= deblank(setstr(fread(fid,8,'char'))');
   d.hdr.unused1			= fread(fid,1,'int16');
   d.hdr.datatype			= fread(fid,1,'int16');
   d.hdr.bitpix				= fread(fid,1,'int16');
   d.hdr.dim_un0			= fread(fid,1,'int16');
   d.hdr.pixdim				= fread(fid,8,'float');
   d.hdr.vox_offset			= fread(fid,1,'float');
   d.hdr.funused1			= fread(fid,1,'float');
   d.hdr.funused2			= fread(fid,1,'float');
   d.hdr.funused3			= fread(fid,1,'float');
   d.hdr.cal_max			= fread(fid,1,'float');
   d.hdr.cal_min			= fread(fid,1,'float');
   d.hdr.compressed			= fread(fid,1,'int32');
   d.hdr.verified			= fread(fid,1,'int32');
   d.hdr.glmax				= fread(fid,1,'int32');
   d.hdr.glmin				= fread(fid,1,'int32');
   
   % read (struct) data_history
   %---------------------------------------------------------------------------
   fseek(fid,148,'bof');
   
   d.hdr.descrip			= deblank(setstr(fread(fid,80,'char'))');
   d.hdr.aux_file			= deblank(setstr(fread(fid,24,'char'))');
   d.hdr.orient				= fread(fid,1,'char');
   d.hdr.origin				= fread(fid,5,'uint16');
   d.hdr.generated			= deblank(setstr(fread(fid,10,'char'))');
   d.hdr.scannum			= deblank(setstr(fread(fid,10,'char'))');
   d.hdr.patient_id			= deblank(setstr(fread(fid,10,'char'))');
   d.hdr.exp_date			= deblank(setstr(fread(fid,10,'char'))');
   d.hdr.exp_time			= deblank(setstr(fread(fid,10,'char'))');
   d.hdr.hist_un0			= deblank(setstr(fread(fid,3,'char'))');
   d.hdr.views				= fread(fid,1,'int32');
   d.hdr.vols_added			= fread(fid,1,'int32');
   d.hdr.start_field		= fread(fid,1,'int32');
   d.hdr.field_skip			= fread(fid,1,'int32');
   d.hdr.omax				= fread(fid,1,'int32');
   d.hdr.omin				= fread(fid,1,'int32');
   d.hdr.smax				= fread(fid,1,'int32');
   d.hdr.smin				= fread(fid,1,'int32');
   
   fclose(fid);
   
   % Put important information in main structure
   %---------------------------------------------------------------------------
   
   d.dim    	  			= d.hdr.dim(2:5)';
   vox 						= d.hdr.pixdim(2:5)';
   if 	vox(4)==0 
      vox(4)=[];
   end
   d.vox       				= vox;
   d.vox_units       		= d.hdr.vox_units;
   d.vox_offset	    		= d.hdr.vox_offset;
   scale     				= d.hdr.funused1;
   d.scale     			  	= ~scale + scale;
   d.global					= [d.hdr.glmin d.hdr.glmax];
   d.calib					= [d.hdr.cal_min d.hdr.cal_max];
   d.calib_units			= d.hdr.cal_units;
   d.origin    				= d.hdr.origin(1:3)';
   d.descrip   				= d.hdr.descrip(1:max(find(d.hdr.descrip)));
   
   if ~isnan(str2double(d.hdr.scannum));    d.df(1)  =str2double(d.hdr.scannum);    end
   if ~isnan(str2double(d.hdr.patient_id)); d.df(2)  =str2double(d.hdr.patient_id); end
   if ~isnan(str2double(d.hdr.exp_date));   d.fwhm(1)=str2double(d.hdr.exp_date);   end
   if ~isnan(str2double(d.hdr.exp_time));   d.fwhm(2)=str2double(d.hdr.exp_time);   end
   if ~isnan(str2double(d.hdr.hist_un0));   d.nconj  =str2double(d.hdr.hist_un0);   end
   
   switch d.hdr.datatype
   case 1
      d.precision 	= 'uint1';
      d.hdr.byte 	= 0;
   case 2
      d.precision 	= 'uint8';
      d.hdr.byte 	= 1;
   case 4
      d.precision 	= 'int16';
      d.hdr.byte 	= 2;
   case 8
      d.precision 	= 'int32';
      d.hdr.byte 	= 4;
   case 16
      d.precision 	= 'float';
      d.hdr.byte 	= 4;
   case 64
      d.precision 	= 'double';
      d.hdr.byte 	= 8;
   otherwise
      d.precision 	= 'Unknown';
      d.hdr.byte 	= 0;
   end
   
else
   d=[];
   errordlg('Check Header File: Existence, Permissions ?','Read Error'); 
end

return
end


%%

function [d]=fmris_read_minc(file,Z,T)

fid=fopen(file);
file=fopen(fid);
fclose(fid);

d.file_name = file;

d.dim 	= fliplr(miinquire(file,'imagesize')');
d.dim 	= d.dim + (d.dim == 0);

if nargin == 1 
   Z = 1:d.dim(3);
   T = 1:d.dim(4);
end

if (T~=0)&(Z~=0)
   if d.dim(4)==1
      images 	= mireadimages(file,Z-1);
   else
      n=length(T);
      images=[];
      for i=0:floor((n-1)/256)
         frames=T((i*256+1):min((i+1)*256,n));
         images 	= [images mireadimages(file,Z-1,frames-1)];
      end
   end
   d.data 		= reshape(images,[d.dim(1:2) length(Z) length(T)]);
   d.calib		= [min(min(min(min(d.data)))) max(max(max(max(d.data))))];
end

[x_vox, y_vox, z_vox] =  miinquire(file,'attvalue','xspace','step','attvalue','yspace','step','attvalue','zspace','step');
d.vox 					= [x_vox y_vox z_vox];
d.vox_units				= '';
d.vox_offset			= 0;
d.precision				= '';
d.calib_units			= '';
[x_origin, y_origin, z_origin] = miinquire(file,'attvalue','xspace','start','attvalue','yspace','start','attvalue','zspace','start');
d.origin 				= [x_origin y_origin z_origin];
d.descrip				= '';

df=miinquire(file, 'attvalue', 'df', 'df');
if ~isempty(df)
   d.df=df;
end

nconj=miinquire(file, 'attvalue', 'nconj', 'nconj');
if ~isempty(nconj)
   d.nconj=nconj;
end

fwhm=miinquire(file, 'attvalue', 'fwhm', 'fwhm');
if ~isempty(fwhm)
   d.fwhm=fwhm;
end

return;
end


%%

function [d]=fmris_read_nifti(file,Z,T)

d.file_name=file;
machineformats='nlbdgcas';
for i=1:length(machineformats)
  fid=fopen(file,'r',machineformats(i));
  if fid == -1;
     fprintf('Invalid file %s',file)
  end
  sizeof_hdr=fread(fid,1,'int');
  if sizeof_hdr==348;
     %Must have found a formt that works, so break and continue using this format
     break;
  else
    %Need to close if file is not native format
    %else if the file format is 's', then 7 stranded files are left orphaned and never closed.
    fclose(fid);
  end
end
data_type=fread(fid,10,'char');
db_name=fread(fid,18,'char');
extents=fread(fid,1,'int');
session_error=fread(fid,1,'short');
regular=char(fread(fid,1,'char')');
dim_info=char(fread(fid,1,'char')');
dim=fread(fid,8,'short');
intent_p =fread(fid,2,'float');
intent_q =fread(fid,2,'uint16');
intent_code =fread(fid,1,'short');
datatype=fread(fid,1,'short');
bitpix=fread(fid,1,'short');
slice_start=fread(fid,1,'short');
pixdim=fread(fid,8,'float');
vox_offset=fread(fid,1,'float');
scl_slope =fread(fid,1,'float');
scl_inter =fread(fid,1,'float');
slice_end=fread(fid,1,'short');
slice_code =char(fread(fid,1,'char')');
xyzt_units =char(fread(fid,1,'char')');
cal_max=fread(fid,1,'float');
cal_min=fread(fid,1,'float');
slice_duration=fread(fid,1,'float');
toffset=fread(fid,1,'float');
glmax=fread(fid,1,'int');
glmin=fread(fid,1,'int');
descrip=char(fread(fid,80,'char')');
aux_file=char(fread(fid,24,'char')');
qform_code =fread(fid,1,'short');
sform_code =fread(fid,1,'short');
quatern_b =fread(fid,1,'float');
quatern_c =fread(fid,1,'float');
quatern_d =fread(fid,1,'float');
qoffset_x =fread(fid,1,'float');
qoffset_y =fread(fid,1,'float');
qoffset_z =fread(fid,1,'float');
srow_x =fread(fid,4,'float');
srow_y =fread(fid,4,'float');
srow_z =fread(fid,4,'float');
intent_name=char(fread(fid,16,'char')');
magic =char(fread(fid,4,'char')');

d.machineformat=machineformats(i);
d.dim=ones(1,4);
if dim(1)==5
   if dim(5)==1 
      dim=[4; dim(2:4); dim(6); zeros(3,1)]
   else
      dim
      fclose(fid);
      return
   end
end
d.dim(1:dim(1))=dim((1:dim(1))+1);
d.vox=zeros(1,3);
d.vox(1:dim(1))=pixdim((1:dim(1))+1);
%Attempt to fill out more information for a complete nifti description
  d.vox_offset  = vox_offset;
  d.scale       = scl_slope;
  d.intercept   = scl_inter;
  d.global      = [glmin glmax];
  d.calib       = [cal_min cal_max];
  if qform_code>0;
    d.origin    = [qoffset_x qoffset_y qoffset_z];
  else
    d.origin    = [srow_x(4) srow_y(4) srow_z(4)];
  end
  d.descrip     = descrip;
d.parent_file=file;

if intent_p(1)>0; d.df(1)=intent_p(1); end;
if intent_p(2)>0; d.df(2)=intent_p(2); end;

intent_q=intent_q/2^16*100;
if intent_q(1)>0; d.fwhm(1)=intent_q(1); end;
if intent_q(2)>0; d.fwhm(2)=intent_q(2); end;

if nargin<2
   Z=1:d.dim(3);
   T=1:d.dim(4);
end

if Z(1)==0 & T(1)==0
   fclose(fid);
   return
end

types=lower(['UINT8  ';'INT16  ';'INT32  ';'FLOAT32';'FLOAT64']);
codes=[2; 4; 8; 16;  64];
d.precision=[];
for i=1:5
   if datatype==codes(i)
      d.precision=deblank(types(i,:));
   end
end
if isempty(d.precision)
   unknown_datatype=datatype
   fclose(fid);
   return
end

%d.byte=bitpix/8
d.byte=4;
if all(Z>0)&all(Z<=d.dim(3))&all(T>0)&all(T<=d.dim(4))
   d.data=zeros(d.dim(1),d.dim(2),length(Z),length(T));
   for t=1:length(T)
      for z=1:length(Z)
         position=d.byte*((T(t)-1)*prod(d.dim(1:3))+(Z(z)-1)*prod(d.dim(1:2)))+vox_offset;
         status=fseek(fid,position,'bof');
         d.data(:,:,z,t)=fread(fid,[d.dim(1) d.dim(2)],d.precision);
      end
   end
  fclose(fid);
else
   Z
   T
   fclose(fid);
   return
end

if scl_slope~=0
   d.data=d.data*scl_slope+scl_inter;
end


return;
end

