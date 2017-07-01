clear all
endian = 'ieee-le';
pth = 'C:\bloodA\Image_analysis_files\examples\SPM analysis\Brain_atlas_Hammers2003';
[nme,pth] = uigetfile(fullfile(pth,'*.img'),'Select dynamic Analyze formatted image file to analyze','Multiselect','off');
img_file = fullfile(pth,nme);
[pth,nme] = fileparts(img_file);
hdr = HDRread(fullfile(pth,nme),endian);
nXY = prod(hdr.dim(2:3));
nZ = hdr.dim(4);
%hdr.datatype = 2;
%hdr.data_type = 2;

for ipl = 1:nZ
Y(ipl,:) = read_plane_static(fullfile(pth,[nme '.img']),hdr,ipl,endian);
end

%regions = [22:23 64:67]; 
%name = 'occipital';

%regions = [28:29 56:59]; 
%name = 'frontal';

%regions = [32:33 60:63]; 
%name = 'parietal';

%regions = [44]; 
%name = 'cc';

%regions = [26 27]; 
%name = 'pcc';

%regions = [58 59]; 
%name = 'mpfc';

%regions = [60 61]; 
%name = 'PcG';

%regions = [34 35]; 
%name = 'caud';

regions = [78 79]; 
name = 'subcall';


Yr = Y .* ismember(Y,regions);

img_write_static(fullfile(pth,[nme '_' name '.img']),permute(reshape(Yr,[hdr.dim(4) hdr.dim(2:3)]),[2 3 1]),hdr,endian,nifti);
