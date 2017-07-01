conntype = 'dwpli';

ftdat = load(['W:\Data\CRPS_resting\EEG\ft_' conntype '\H07.100Hz.Exp3_BELB_subcomp_ft' conntype '.mat'],'matrix');
eldat = load(['W:\Data\CRPS_resting\EEG\eeglab_' conntype '\H07.100Hz.Exp3_BELB_subcomp_eeglab' conntype '.mat'],'matrix');
ftdat  = ftdat.matrix;
eldat  = eldat.matrix;


ftdat1 = squeeze(ftdat(1,:,:))
eldat1 = squeeze(eldat(1,:,:))

figure
for f = 1:5

    subplot(5,2,(f-1)*2+1)
    imagesc(squeeze(ftdat(f,:,:)))
    colormap(jet)
    axis image
    title(['ft.' conntype '.freq' num2str(f)])

    subplot(5,2,(f-1)*2+2)
    imagesc(squeeze(eldat(f,:,:)))
    colormap(jet)
    axis image
    title(['eeglab.' conntype '.freq' num2str(f)])

end
