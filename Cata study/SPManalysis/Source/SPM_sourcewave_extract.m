function [Ds Sw] = SPM_sourcewave_extract(fname, val, Cnii, sep)

D = spm_eeg_load(fname);

% get MNI corrdinates of Cnii cluster
img=spm_read_vols(Cnii);
ind3={};
if sep
    uimg = unique(img)';
    uimg(uimg==0) = [];
    for u = uimg
        c = find(img==u);
        [x y z] = ind2sub(size(img),c);
        ind3{u} = [x y z];
    end
else
    c = find(img>0);
    [x y z] = ind2sub(size(img),c);
    ind3{1} = [x y z];
end

for i = 1:length(ind3)
    disp(['sub-cluster ' num2str(i) '/' num2str(length(ind3))])
    mni{i} = cor2mni(ind3{i}, Cnii.mat); % convert to MNI
    D.val=val;
    D.inv{val}.source.XYZ  = mni{i}; %- (n x 3) matrix of MNI coordinates
    D.inv{val}.source.cluster = 1; % combine all coords into a single source cluster
    D.inv{val}.source.rad = 0;   %- radius (mm) of VOIs (default 5 mm)
    %D.inv{val}.source.label %- label(s) for sources (cell array)
    D.inv{val}.source.type = 'evoked';  %- output type ('evoked'/'trials')
    D.inv{val}.source.save = 0;
    D.inv{val}.source.fname = ''; %- output file name

    Ds = spm_eeg_inv_extract_CAB(D);
    Sw{i} = Ds(:,:,:);
end

% plot first condition as an example
%figure;plot(Ds(:,:,1))
