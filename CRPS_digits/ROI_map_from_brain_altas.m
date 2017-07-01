% NEEDS RE-WRITING USING IMCALC
clear all
endian = 'ieee-le';
pth = 'W:\Brain_atlas_Hammers2003';
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

regions = {'L_postcen',60;
	'R_postcen',61;
	'L_sup_par',62;
	'R_sup_par',63;
	'L_inf_par',32;
	'R_inf_par',33;
    };
Nreg = size(regions,1);

for r = 1:Nreg
    region = regions{r,2};
    name = regions{r,1};

    Yr = Y .* ismember(Y,region);

    img_write_static(fullfile(pth,[nme '_' name '.img']),permute(reshape(Yr,[hdr.dim(4) hdr.dim(2:3)]),[2 3 1]),hdr,endian,nifti);
    
end
