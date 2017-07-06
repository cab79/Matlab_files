% 1. loads a topography*time image
% 2. finds extreme (positive or negative, not both) values within a time
% window of interest; defined by the number of standard deviations from the
% mean
% 3. creates Maximum Intensity Projection (collapsing topography over the 
% time window of interest) to use as a mask.

clear all

% path to image files
S.path = 'C:\Users\Chris\Google Drive\SPM\SPM Scripts Sarah Version\t-4000_500_b-3750_-3500_mspm12_cS1HC11_Session1_AppliedICA_2ndRej';

% image file name(s) including .nii extension - can be move than one
S.fnames = {
    'condition_S4.nii'
    };

% time limits of input image file (ms)
S.timelim = [-4000 500];

% sampling rate of input image file (Hz)
S.sr = 500;

% time window to mask (inclusive) - all other latencies will be zero
S.timewin = [-500 0];

% number of standard deviations - e.g. 2 (positive values more than 2 x SDs) 
% or -2 (negative values less than -2 x SDs) from the mean
% If the mask area is not large enough, make this value smaller
S.nSD = 1;

%% RUN

% load SPM EEG data file to find latency info and sample rate
lats=S.timelim(1):1000/S.sr:S.timelim(2);

nfiles = length(S.fnames);
for f = 1:nfiles
    % load image
    nii = load_nii(fullfile(S.path,S.fnames{f}));
    % check it has the correct number of time points
    dim_img=size(nii.img);
    if dim_img(3)~=length(lats)
        dbstop if error
        error('number of data points in nii image do not match S.timelim and S.sr')
    end
    % find latency index of S.timewin
    lat_ind = find(lats==S.timewin(1)):find(lats==S.timewin(2));
    % extract data from timewin
    img_tw=nii.img(:,:,lat_ind);
    % create mean image
    meanimg = nanmean(img_tw,3); % mean over time
    mmimg = nanmean(meanimg(:)); % mean of mean img
    sdimg = nanstd(meanimg(:)); % sd of mean img
    % threshold mean image to S.nSD
    if S.nSD>0
        img_thresh = double(meanimg>(mmimg+S.nSD*sdimg));
    elseif S.nSD<0
        img_thresh = double(meanimg<(mmimg+S.nSD*sdimg));
    end
    % create mask
    mask_img = zeros(dim_img);
    mask_img(:,:,lat_ind) = repmat(img_thresh,1,1,length(lat_ind));
    nii.img=mask_img;
    % save image
    [~,nme,~]=fileparts(S.fnames{f});
    save_nii(nii,fullfile(S.path,[nme '_thresh_SD' num2str(S.nSD) '_t' num2str(S.timewin(1)) '_t' num2str(S.timewin(2)) '.nii']));
end

