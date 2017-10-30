function [Ds Sw] = SPM_sourcewave_extract(fname, val, Cnii, sep, voi, use_VOI)

D = spm_eeg_load(fname);

% get MNI corrdinates of Cnii cluster

if sep
    img=spm_read_vols(Cnii);
    ind3={};
    uimg = unique(img)';
    uimg(uimg==0) = [];
    for u = uimg
        c = find(img==u);
        [x y z] = ind2sub(size(img),c);
        ind3{u} = [x y z];
        mni{u} = cor2mni(ind3{u}, Cnii.mat); % convert to MNI
    end
elseif use_VOI
    mni = {voi.xY.XYZmm'};
else
    c = find(img>0);
    [x y z] = ind2sub(size(img),c);
    ind3{1} = [x y z];
    mni{1} = cor2mni(ind3{1}, Cnii.mat); % convert to MNI
end

for i = 1:length(mni)
    disp(['sub-cluster ' num2str(i) '/' num2str(length(mni))])
    D.val=val;
    D.inv{val}.source.XYZ  = mni{i}; %- (n x 3) matrix of MNI coordinates
    D.inv{val}.source.cluster = 1; % combine all coords into a single source cluster
    D.inv{val}.source.rad = 0;   %- radius (mm) of VOIs (default 5 mm)
    %D.inv{val}.source.label %- label(s) for sources (cell array)
    D.inv{val}.source.type = 'evoked';  %- output type ('evoked'/'trials')
    D.inv{val}.source.save = 0;
    D.inv{val}.source.fname = ['temp' num2str(val)]; %- output file name

    Ds = spm_eeg_inv_extract_CAB(D);
    Sw{i} = Ds(:,:,:);
end

% plot first condition as an example
%figure;plot(Ds(:,:,1))
