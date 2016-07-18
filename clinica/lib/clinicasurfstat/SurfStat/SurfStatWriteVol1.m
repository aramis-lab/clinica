function SurfStatWriteVol( d, Z, T );

%Writes volumetric image data in MINC, ANALYZE, NIFTI or AFNI format. 
% 
% Usage: fmris_write_image( d [, Z, T] ). 
%
% d.file_name = file name with extension .mnc, .img, .nii or .brik as above
% Z           = vector of slices.
% T           = vector of times. 
% If Z and T are omitted, writes the entire volume.

file = deblank(d.file_name);

[base,ext]=fileparts2(file);
d.file_name=[base ext];
   
switch lower(ext)
case '.mnc'
   if nargin == 3
      if length(T)<=160
         fmris_write_minc(d,Z,T);
      else
         fn=deblank(d.file_name);
         d.file_name=[fn(1:(length(fn)-3)) 'nii'];
         fmris_write_nifti(d,Z,T);
      end
   else
      if d.dim(4)<=160
         fmris_write_minc(d);
      else
         fn=deblank(d.file_name);
         d.file_name=[fn(1:(length(fn)-3)) 'nii'];
         fmris_write_nifti(d);
      end
   end
case '.img'
   if nargin == 3
      fmris_write_analyze(d,Z,T);
   else
      fmris_write_analyze(d);
   end
case '.brik'
   if nargin == 3
      fmris_write_afni(d,Z,T);
   else
      fmris_write_afni(d);
   end
case '.nii'
   if nargin == 3
      fmris_write_nifti(d,Z,T);
   else
      fmris_write_nifti(d);
   end
otherwise
   ['Unknown file extension']
end

return
end

%%

function [base,ext]=fileparts2(string)
if isstr(string)
    [path,name,ext]=fileparts(deblank(string));
else
    [path,name,ext]=fileparts(deblank(string.file_name));
end
if strcmp(ext,'.gz')
   [path2,name,ext]=fileparts(name);
end   
if isempty(path)
   base=name;
else
   base=[path '/' name];
end

return
end


%%

function	[d]=fmris_write_afni(d,Z,T);

existinfo=0;
if ~isfield(d,'dim') & isfield(d,'parent_file')
   [path,name,ext]=fileparts(deblank(d.parent_file));
   [err, Infoparent] = BrikInfo([path '/' name '.HEAD']);
   d.dim=[Infoparent.DATASET_DIMENSIONS(1:3) Infoparent.DATASET_RANK(2)];
   existinfo=1;
end

if nargin==1 
   Z=1:d.dim(3);
   T=1:d.dim(4);
end

Opt.Prefix=d.file_name(1:(length(d.file_name)-5));
Opt.Slices=Z;
Opt.Frames=T;
Opt.NoCheck=1;

if Z(1)==1 & T(1)==1
   if ~existinfo
      [path,name,ext]=fileparts(deblank(d.parent_file));
      [err, Infoparent] = BrikInfo([path '/' name '.HEAD']);
   end
   [path,name,ext]=fileparts(deblank(d.file_name));
   Info.LABEL_1=name;
   Info.DATASET_NAME=['./' name];
   if isfield(d,'origin')
      Info.ORIGIN=d.origin;
   else
      Info.ORIGIN=Infoparent.ORIGIN;
   end
   if isfield(d,'vox')
      Info.DELTA=d.vox;
   else
      Info.DELTA=Infoparent.DELTA;
   end
   Info.SCENE_DATA=Infoparent.SCENE_DATA;
   Info.ORIENT_SPECIFIC=Infoparent.ORIENT_SPECIFIC;
   Info.TYPESTRING=Infoparent.TYPESTRING;
   Opt.NoCheck=0;
end

Info.DATASET_DIMENSIONS=[d.dim(1:3) 0 0];
Info.DATASET_RANK=[3 d.dim(4) 0 0 0 0 0 0];
Info.BRICK_TYPES=repmat(3,1,d.dim(4));
Info.TypeName='float';
Info.TypeBytes=4;
Info.BYTEORDER_STRING='MSB_FIRST';
Info.MachineFormat='ieee-be';

if isfield(d,'df')
   if ~isempty(d.df)
      Info.WORSLEY_DF=d.df;
   end
end

if isfield(d,'nconj')
   if ~isempty(d.nconj)
      Info.WORSLEY_NCONJ=d.nconj;
   end
end

if isfield(d,'fwhm')
   if ~isempty(d.fwhm)
      Info.WORSLEY_FWHM=d.fwhm;
   end
end

[err, ErrMessage, Info] = WriteBrik(d.data, Info, Opt);

return
end


%%

function [d]=fmris_write_analyze(d,Z,T);

% Write Analyze format Image Volumes
%
% Usage
%
% rpm_write_analyze(d,Z,T);
% 
% Writes image data and attributes from the structure d.
% The following fields are required from d. 
%
% d.file_path: '/kop1/data/fdg/'
% d.file_name: 'n03309_3d_dy2_CS-Cal_K1'
% d.data: [128x128x31 double]
% d.vox: [2.0900 2.0900 3.4200]
% d.vox_units: 'mm'
% d.vox_offset: 0
% d.calib_units: 'min^-1'
% d.origin: [0 0 0];
% d.descrip: 'Parametric Image - K1 (1)'
%
% All other information is generated automatically.
%
% (c) Roger Gunn & John Aston  

if isfield(d,'parent_file')
   d2=SurfStatReadVol(d.parent_file,0,0);
   d.vox=d2.vox;
   d.vox_units=d2.vox_units;
   d.calib_units='';
   d.origin=d2.origin;
   d.vox_offset=d2.vox_offset;
end
if ~isfield(d,'descrip')
   d.descrip='';
end
file=d.file_name;

d.origin=-d.origin./d.vox(1:3);

if length(size(d.data))<5&size(d.data,1)>1&size(d.data,2)>1&length(d.origin)<4
    
    d.calib = [1 1];
    d.precision='float';  
    
    % Write Header
    if nargin==1
        [d]=Write_Analyze_Hdr(d);
    elseif nargin==3
        if (T(1)==1&Z(1)==1)
            [d]=Write_Analyze_Hdr(d);
        else
            d3 = fmris_read_image(file,0,0);
            d.hdr = d3.hdr;
        end      
    end
    
    if ~isstruct(d);
        return
    end
    % Write Image 
    if nargin==1
        if ~isempty(d)
            fid = fopen(file,'w','n');
            if fid > -1
                for t=1:d.hdr.dim(5)
                    for z=1:d.hdr.dim(4)
                        fwrite(fid,d.data(:,:,z,t)/d.hdr.funused1,d.precision);
                    end
                end
                fclose (fid);
            else
                errordlg('Cannot open file for writing  ','Write Error');d=[];return;
            end
        end
    elseif nargin==3    
        if T(1)~=1|Z(1)~=1
            if ~exist(file,'file')
                errordlg('Please write Plane 1 Frame 1 first','Write Error');return;
            end
        else
            fid=fopen(file,'w','n');fclose(fid);
        end
        fid = fopen(file,'r+','n');
        if fid > -1
            if T(1)==1&Z(1)==1
                plane=zeros(d.dim(1:2));
                for t=1:d.hdr.dim(5)
                    for z=1:d.hdr.dim(4)
                        fwrite(fid,plane,d.precision);
                    end
                end      
            end
            
            for t=1:length(T)
                for z=1:length(Z)
                    fseek(fid,4*d.dim(1)*d.dim(2)*((T(t)-1)*d.dim(3)+Z(z)-1),'bof');
                    if length(Z)~=1;
                        fwrite(fid,d.data(:,:,z,t),'float');
                    else
                        fwrite(fid,d.data(:,:,t),'float');
                    end
                end
            end
            
            fclose (fid);
            
        else
            errordlg('Cannot open file for writing  ','Write Error');d=[];return;
            
        end
    end
    
else
    errordlg('Incompatible data structure: Check dimension and Origin  ','Write Error'); 
end

return;
end

%%

function [d]=Write_Analyze_Hdr(d);

% Write Analyze Header from the structure d
% Adapted from John Ashburners spm_hwrite.m

d.file_name_hdr=[d.file_name(1:(length(d.file_name)-3)) 'hdr'];
file=d.file_name_hdr;

fid   			= fopen(file,'w','n');
if fid > -1
    d.hdr.data_type 			= ['dsr      ' 0];
    d.hdr.db_name	  		= ['                 ' 0];
    if isfield(d,'dim')
        d.hdr.dim    			= [4 1 1 1 1 0 0 0];
        d.hdr.dim(2:(1+length(d.dim(find(d.dim)))))= d.dim(find(d.dim));
    else
        d.hdr.dim    			= [4 1 1 1 1 0 0 0];
        d.hdr.dim(2:(1+length(size(d.data))))   = size(d.data);
    end   
    
    d.hdr.pixdim 			= [4 0 0 0 0 0 0 0];
    d.hdr.pixdim(2:(1+length(d.vox))) = d.vox;
    d.hdr.vox_units			= [0 0 0 0];
    d.hdr.vox_units(1:min([3 length(d.vox_units)])) = d.vox_units(1:min([3 length(d.vox_units)]));
    d.hdr.vox_offset 		= d.vox_offset;
    d.hdr.calmin				= d.calib(1);
    d.hdr.calmax				= d.calib(2);
    switch d.precision
    case 'uint1'  % 1  bit
        d.hdr.datatype 		= 1;
        d.hdr.bitpix 			= 1;
        d.hdr.glmin			= 0;
        d.hdr.glmax 			= 1;
        d.hdr.funused1		= 1;   
    case 'uint8'  % 8  bit
        % d.hdr.datatype 		= 2;
        % d.hdr.bitpix 			= 8;
        % d.hdr.glmin 			= 0;
        % d.hdr.glmax 			= 255;
        % d.hdr.funused1	= abs(d.hdr.calmin)/255;
        errordlg('You should write a float image','8 Bit Write Error');d=[];return;
    case 'int16'  % 16 bit
        d.hdr.datatype 		= 4;
        d.hdr.bitpix  		= 16;
        if abs(d.hdr.calmin)>abs(d.hdr.calmax)
            d.hdr.funused1  	= abs(d.hdr.calmin)/(2^15-1);
        else
            d.hdr.funused1	= abs(d.hdr.calmax)/(2^15-1);
        end
        d.hdr.glmin 			= round(d.hdr.funused1*d.hdr.calmin);
        d.hdr.glmax 			= round(d.hdr.funused1*d.hdr.calmin);
    case 'int32'  % 32 bit
        d.hdr.datatype 		= 8;
        d.hdr.bitpix  		= 32;
        if abs(d.hdr.calmin)>abs(d.hdr.calmax)
            d.hdr.funused1  	= abs(d.hdr.calmin)/(2^31-1);
        else
            d.hdr.funused1	= abs(d.hdr.calmax)/(2^31-1);
        end
        d.hdr.glmin 			= round(d.hdr.funused1*d.hdr.calmin);
        d.hdr.glmax 			= round(d.hdr.funused1*d.hdr.calmin);
    case 'float'  % float  (32 bit)
        d.hdr.datatype 		= 16;
        d.hdr.bitpix 	 		= 32;
        d.hdr.glmin 			= 0;
        d.hdr.glmax 			= 0;
        d.hdr.funused1 		= 1;
    case 'double' % double (64 bit) 
        d.hdr.datatype 		= 64;
        d.hdr.bitpix  		= 64;
        d.hdr.glmin 			= 0;
        d.hdr.glmax 			= 0;
        d.hdr.funused1 		= 1;
    otherwise
        errordlg('Unrecognised precision (d.type)','Write Error');d=[];return;
    end
    d.hdr.descrip 			= zeros(1,80);d.hdr.descrip(1:min([length(d.descrip) 79]))=d.descrip(1:min([length(d.descrip) 79]));
    d.hdr.aux_file        	= ['none                   ' 0];
    d.hdr.origin          	= [0 0 0 0 0];d.hdr.origin(1:length(d.origin))=d.origin;
    
    
    % write (struct) header_key
    %---------------------------------------------------------------------------
    fseek(fid,0,'bof');
    
    fwrite(fid,348,					'int32');
    fwrite(fid,d.hdr.data_type,	'char' );
    fwrite(fid,d.hdr.db_name,		'char' );
    fwrite(fid,0,					'int32');
    fwrite(fid,0,					'int16');
    fwrite(fid,'r',					'char' );
    fwrite(fid,'0',					'char' );
    
    
    % write (struct) image_dimension
    %---------------------------------------------------------------------------
    fseek(fid,40,'bof');
    
    fwrite(fid,d.hdr.dim,			'int16');
    fwrite(fid,d.hdr.vox_units,	'char' );
    fwrite(fid,zeros(1,8),			'char' );
    fwrite(fid,0,					'int16');
    fwrite(fid,d.hdr.datatype,	'int16');
    fwrite(fid,d.hdr.bitpix,		'int16');
    fwrite(fid,0,					'int16');
    fwrite(fid,d.hdr.pixdim,		'float');
    fwrite(fid,d.hdr.vox_offset,	'float');
    fwrite(fid,d.hdr.funused1,	'float');
    fwrite(fid,0,					'float');
    fwrite(fid,0,					'float');
    fwrite(fid,d.hdr.calmax,		'float');
    fwrite(fid,d.hdr.calmin,		'float');
    fwrite(fid,0,					'int32');
    fwrite(fid,0,					'int32');
    fwrite(fid,d.hdr.glmax,		'int32');
    fwrite(fid,d.hdr.glmin,		'int32');
    
    % write (struct) data_history
    %---------------------------------------------------------------------------
    fwrite(fid,d.hdr.descrip,		'char');
    fwrite(fid,d.hdr.aux_file,   	'char');
    fwrite(fid,0,           		'char');
    fwrite(fid,d.hdr.origin,     'uint16');
    
    fwrite(fid,zeros(1,10), 		'char');
    
    if isfield(d,'df')
       fwrite(fid,sprintf('%10g',d.df(1)),'char');
       if length(d.df)>1
          fwrite(fid,sprintf('%10g',d.df(2)),'char');
       else
          fwrite(fid,zeros(1,10),'char');
       end
    else
       fwrite(fid,zeros(1,20),'char');
    end
    
    if isfield(d,'fwhm')
       fwrite(fid,sprintf('%10g',d.fwhm(1)),'char');
       if length(d.fwhm)>1
          fwrite(fid,sprintf('%10g',d.fwhm(2)),'char');
       else
          fwrite(fid,zeros(1,10),'char');
       end
    else
       fwrite(fid,zeros(1,20),'char');
    end
    
    if isfield(d,'nconj')
       fwrite(fid,sprintf('%3g',d.nconj),'char');
    else
       fwrite(fid,zeros(1,3),'char');
    end
    
    fwrite(fid,zeros(1,32),'char');
    
    s   = ftell(fid);
    fclose(fid);
else
    errordlg('Cannot open file for writing  ','Write Error');d=[];return;
end

return;
end


%%

function	[d]=fmris_write_minc(d,Z,T);

if ~(d.file_name(1)=='/' | d.file_name(2)==':')
   if d.file_name(1:2)=='./'
      d.file_name=[pwd d.file_name(3:length(file))];
   else
      d.file_name=[pwd '/' d.file_name];
   end
end

fid=fopen(d.parent_file);
d.parent_file=fopen(fid);
fclose(fid);

if ~isfield(d,'dim') 
   d.dim 	= fliplr(miinquire(d.parent_file,'imagesize')');
   d.dim 	= d.dim + (d.dim == 0);
end

if nargin==1
   Z=1:d.dim(3);
   T=1:d.dim(4);
end

if isfield(d,'precision')
   precision=d.precision;
   if isempty(precision)
       precision='float';
   end
else
   precision='float';
end

if T(1)==1 & Z(1)==1
   dim=d.dim;
   if dim(4)==1
      dim(4)=0;
   end
   if exist(d.file_name,'file')
      delete(d.file_name);
   end
   newh=newimage(d.file_name,fliplr(dim),d.parent_file,precision);
   if isfield(d,'df')
      miwriteatt(d.file_name,'df','df',d.df);
   end
   if isfield(d,'nconj')
      miwriteatt(d.file_name,'nconj','nconj',d.nconj);
   end
   if isfield(d,'fwhm')
      miwriteatt(d.file_name,'fwhm','fwhm',d.fwhm);
   end
else
   newh=openimage(d.file_name,'w');
end

d.data=squeeze(reshape(d.data,d.dim(1)*d.dim(2),length(Z),length(T)));

if d.dim(4)<=1
   putimages(newh,d.data,Z);
elseif length(T)==1|length(Z)==1
   putimages(newh,d.data,Z,T);
else
   for i=1:length(T)
      putimages(newh,squeeze(d.data(:,:,i)),Z,T(i));
   end
end

closeimage(newh);

return
end


%%

function [d]=fmris_write_nifti(d,Z,T)

if nargin<2
   Z=1;
   T=1;
end

if Z(1)==1 & T(1)==1
   fidout=fopen(d.file_name,'w');
   isniftiparent=0;
   if isfield(d,'parent_file')
      if exist(d.parent_file)
         pf=deblank(d.parent_file);
         ext=pf((length(pf)-2):length(pf));
         isniftiparent=(ext=='nii');
      end
   end
   if isniftiparent
      machineformats='nlbdgcas';
      for i=1:length(machineformats)
         fid=fopen(d.parent_file,'r',machineformats(i));
         sizeof_hdr=fread(fid,1,'int');
         if sizeof_hdr==348;
            %Must have found a formt that works, so break and continue using this format
            break
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
      intent_p =fread(fid,3,'float');
      intent_code =fread(fid,1,'short');
      datatype=fread(fid,1,'short');
      bitpix=fread(fid,1,'short');
      slice_start=fread(fid,1,'short');
      pixdim=fread(fid,8,'float');
      vox_offset=fread(fid,1,'float');
      scl_slope=fread(fid,1,'float');
      scl_inter=fread(fid,1,'float');
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
      magic =fread(fid,4,'char');
      fclose(fid);
   else
      data_type='          ';
      db_name='                  ';
      extents=0;
      session_error=0;
      regular='r';
      dim_info=' ';
      slice_start=0;
      pixdim=ones(1,8);
      slice_code =' ';
      xyzt_units =' ';
      cal_max=25500;
      cal_min=3;
      slice_duration=0;
      toffset=0;
      aux_file='                        ';
      qform_code =1;
      sform_code =0;
      quatern_b =0;
      quatern_c =1;
      quatern_d =0;
      qoffset_x =d.origin(1);
      qoffset_y =d.origin(2);
      qoffset_z =d.origin(3);
      srow_x =[0 0 0 0];
      srow_y =[0 0 0 0];
      srow_z =[0 0 0 0];
      intent_name='                ';
   end
   
   sizeof_hdr=348;
   datatype=16;
   bitpix=32;
   vox_offset=352;
   scl_slope =0;
   scl_inter =0;
   slice_end=0;
   glmax=0;
   glmin=0;
   descrip=['FMRISTAT' repmat(' ',1,72)];
   magic =[double('n+1') 0]';
   
   dim=ones(1,8);      
   dim(1)=max(find(d.dim>1));
   dim((1:dim(1))+1)=d.dim(1:dim(1));
   if isfield(d,'vox')
      pixdim=ones(1,8);
      pixdim(1)=-1;
      pixdim((1:dim(1))+1)=d.vox(1:dim(1));
   end
   
   intent_p=zeros(1,2);
   if isfield(d,'df') 
      intent_p(1:length(d.df))=d.df; 
      intent_code=length(d.df)+2;
   else
      intent_code=0;
   end
   
   intent_q=zeros(1,2);
   if isfield(d,'fwhm'); 
      intent_q(1:length(d.fwhm))=d.fwhm; 
   end;
   intent_q=round(intent_q/100*2^16);
   intent_q=(intent_q.*(intent_q>=0)-2^16).*(intent_q<2^16)+2^16;

   descrip=['FMRISTAT' repmat(' ',1,72)];
   
   fwrite(fidout,sizeof_hdr,'int');
   fwrite(fidout,data_type,'char');
   fwrite(fidout,db_name,'char');
   fwrite(fidout,extents,'int');
   fwrite(fidout,session_error,'short');
   fwrite(fidout,regular,'char');
   fwrite(fidout,dim_info,'char');
   fwrite(fidout,dim,'short');
   fwrite(fidout,intent_p,'float');
   fwrite(fidout,intent_q,'short');
   fwrite(fidout,intent_code,'short');
   fwrite(fidout,datatype,'short');
   fwrite(fidout,bitpix,'short');
   fwrite(fidout,slice_start,'short');
   fwrite(fidout,pixdim,'float');
   fwrite(fidout,vox_offset,'float');
   fwrite(fidout,scl_slope ,'float');
   fwrite(fidout,scl_inter ,'float');
   fwrite(fidout,slice_end,'short');
   fwrite(fidout,slice_code,'char');
   fwrite(fidout,xyzt_units,'char');
   fwrite(fidout,cal_max,'float');
   fwrite(fidout,cal_min,'float');
   fwrite(fidout,slice_duration,'float');
   fwrite(fidout,toffset,'float');
   fwrite(fidout,glmax,'int');
   fwrite(fidout,glmin,'int');
   fwrite(fidout,descrip,'char');
   fwrite(fidout,aux_file,'char');
   fwrite(fidout,qform_code,'short');
   fwrite(fidout,sform_code,'short');
   fwrite(fidout,quatern_b,'float');
   fwrite(fidout,quatern_c,'float');
   fwrite(fidout,quatern_d,'float');
   fwrite(fidout,qoffset_x,'float');
   fwrite(fidout,qoffset_y,'float');
   fwrite(fidout,qoffset_z,'float');
   fwrite(fidout,srow_x,'float');
   fwrite(fidout,srow_y,'float');
   fwrite(fidout,srow_z,'float');
   fwrite(fidout,intent_name,'char');
   fwrite(fidout,magic,'char');
   
   fwrite(fidout,0,'float');
else
   fidout=fopen(d.file_name,'r+');
end

if nargin<2
   Z=1:d.dim(3);
   T=1:d.dim(4);
end

vox_offset=352;
if ~isfield(d,'precision') | isempty(d.precision); d.precision='float32'; end;
if ~isfield(d,'byte') | isempty(d.byte); d.byte=4; end;

if Z(1)==1 & T(1)==1 & ~(all(Z==1:d.dim(3)) & all(T==1:d.dim(4)))
   for t=1:d.dim(4)
      fwrite(fidout,zeros(1,prod(d.dim(1:3))),d.precision);
   end
end

if all(Z>0)&all(Z<=d.dim(3))&all(T>0)&all(T<=d.dim(4))
   for t=1:length(T)
      for z=1:length(Z)
         position=d.byte*((T(t)-1)*prod(d.dim(1:3))+(Z(z)-1)*prod(d.dim(1:2)))+vox_offset;
         status=fseek(fidout,position,'bof');
         if length(size(d.data))==4
            fwrite(fidout,d.data(:,:,z,t),d.precision);
         elseif length(T)==1 
            fwrite(fidout,d.data(:,:,z),d.precision);
         elseif length(Z)==1  
            fwrite(fidout,d.data(:,:,t),d.precision);
         else
            'Slices and/or frames do not match data'
         end
      end
   end
else
   'Slices and/or frames out of range:'
   Z
   T
end

fclose(fidout);

return;
end

