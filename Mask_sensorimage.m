% 1. loads a topography*time image
% 2. finds extreme (positive or negative, not both) values within a time
% window of interest; defined by the number of standard deviations from the
% mean
% 3. creates Maximum Intensity Projection (collapsing topography over the 
% time window of interest) to use as a mask.

clear all

% path to image files
S.path = '';

% image file name(s) including .nii extension - can be move than one
S.fnames = {
    ''
    };

% path and filename of example SPM EEG data file for this study used to generate images
S.spmfile = 'mspm_';

% time window to mask (inclusive) - all other latencies will be zero
S.timewin = [-500 0];

% number of standard deviations - e.g. 2 (positive values more than 2 x SDs) 
% or -2 (negative values less than -2 x SDs) from the mean
S.nSD = 2;

%% RUN

% load SPM EEG data file to find latency info and sample rate
load(S.spmfile)
if ~isstruct(D)
    D=struct(D);
end
lats=;
sr=;
Ntp_data=lats/sr;

nfiles = length(S.fnames);
for f = 1:nfiles
    % load image
    nii = load_nii(fullfile(S.path,S.fnames{f}));
    % check it has the correct number of time points
    dim_img=size(nii.img);
    if dim_img(3)~=Ntp_data
        dbstop if error
        error('number of data points in SPM data file and nii image do not match')
    end
    % find latency index of S.timewin
    lat_ind = find(lats==S.timewin(1)):find(lats==S.timewin(2));
    % extract data from timewin
    img_tw=nii.img(:,:,lat_ind);
    % create mean image
    meanimg = mean(img_tw,3); % mean over time
    mmimg = mean(meanimg(:)); % mean of mean img
    sdimg = std(meanimg(:)); % sd of mean img
    % threshold mean image to S.nSD
    if S.nSD>0
        img_thresh = meanimg(meanimg>mmimg+S.nSD*sdimg);
    elseif S.nSD<0
        img_thresh = meanimg(meanimg<mmimg-S.nSD*sdimg);
    end
    % create mask
    mask_img = zeros(dim_img);
    mask_img(:,:,lat_ind) = repmat(img_thresh,length(lat_ind),1);
    % save image
    [~,nme,~]=fileparts(S.fnames{f});
    savenii(fullfile(S.path,[nme]))
end

