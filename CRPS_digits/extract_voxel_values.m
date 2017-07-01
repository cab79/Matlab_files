% Extracts values from an image data file based on regions/clusters of
% voxels in mask marsbar or image file(s) or voxel.
%
% Usage: [Ym R info] = extract_voxel_values(mask,data)
%
% Ym    - matrix of output mean of values across voxels per image (row) per
%         region (column).
% R     - structure containing voxel info per region (R), per image (I):
%         R.I.Ya  - output vector of values of all voxels.
%         R.I.xyz - matrix indices of voxels in mask region/cluster.
%         R.I.mni - mni coordinates of voxels in mask region/cluster.
%         E.g. to get vector of value for region r, image i, after 
%         extraction, type R(r).I(i).Ya
% info  - structure containing region and image files used:
%         info.regions - string array of region files used.
%         info.images - string array of image files used.
% mask  - string array of mask image file name(s).
% data  - string array of image data file name(s). Can be any analyze or
%         nifti image that contains voxel values of interest, which can be
%         beta values, contrast values, t values etc.
%         e.g. ['beta_0005.img';'beta_0001.img'];
%
% If mask or data not defined, GUI asks to select files.
%
% Example usage:
% extract_voxel_values('L_Fusiform_roi.mat','beta_0005.img')
%
% Dependencies: Marsbar (if *_roi.mat ROI mask files used).
% 
% Modified from extract_voxel_values_old.m by Josh Goh 24 May 2013.
 
function [Ym R info] = extract_voxel_values(mask,data)
 
% Check input
if nargin<1
    mask = spm_select(Inf,'any','Select mask ROI files',[],pwd);
    data = spm_select(Inf,'image','Select data file (*.img or *.nii)',[],pwd);
end
 
info.regions = mask;
info.images  = data;
 
% Loop image
for i = 1:size(data,1)
    
    % Read image header
    V = spm_vol(deblank(data(i,:)));
    
    % Loop regions
    for r = 1:size(mask,1)
        
        % Get voxels from mask
        [~,~,e] = fileparts(deblank(mask(r,:)));
        switch e
            case {'.mat'}
                roi = maroi(deblank(mask(r,:)));
                xyz = voxpts(roi,deblank(data(i,:)));
            case {'.nii','.img'}
                maskdata = spm_read_vols(spm_vol(deblank(mask(r,:))));
                [x,y,z] = ind2sub(size(maskdata),find(maskdata));
                xyz = [x y z]';
        end
        
        % Extract data
        Ya = spm_sample_vol(V,xyz(1,:),xyz(2,:),xyz(3,:),0);
        R(r).I(i).Ya = Ya;
        
        % Compute mean across voxels
        Ym(i,r) = mean(Ya(~isnan(Ya)));
        
        % Compute MNI coordinates
        R(r).I(i).mni = vox2mni(V.mat,xyz);
        R(r).I(i).xyz = xyz;
        
    end
end