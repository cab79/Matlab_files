function EEG = FTrejman(EEG,filterset)
%convert to FT and manually remove chan and trial condidering a particular freq

if isempty(filterset)
    filterset = [0 0];
end

EEGf=EEG;
if filterset(1)>0; EEGf = pop_eegfiltnew( EEG, filterset(1), 0, [], 0);end
if filterset(2)>0; EEGf = pop_eegfiltnew( EEGf, 0, filterset(2), [], 0);end

FT = convertoft(EEGf);
orig_chans=FT.label';
cfg =[];
cfg.method = 'summary';
cfg.alim = 1e-5;
cfg.keepchannel='no';
cfg.keeptrial='yes';
FT = ft_rejectvisual(cfg, FT);
[~, rejchan] = setdiff(orig_chans, FT.label);
[~, ~, rejtrial] = intersect(FT.cfg.artfctdef.summary.artifact(:,1),FT.sampleinfo(:,1));
EEG = eeg_interp(EEG, rejchan);
EEG = pop_select(EEG, 'notrial', rejtrial);