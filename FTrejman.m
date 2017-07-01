function EEG = FTrejman(EEG,filterset,varargin)
%convert to FT and manually remove chan and trial condidering a particular freq

if ~isempty(varargin)
    format = varargin;
else
    format = 'EEGLAB';
end

if isempty(filterset)
    filterset = [0 0];
end

if strcmp(format,'EEGLAB')
    EEGf=EEG;
    if filterset(1)>0; EEGf = pop_eegfiltnew( EEG, filterset(1), 0, [], 0);end
    if filterset(2)>0; EEGf = pop_eegfiltnew( EEGf, 0, filterset(2), [], 0);end
    FT = convertoft(EEGf);
elseif strcmp(format,'SPM')
    FT = EEG.fttimelock;
end
orig_chans=FT.label';
cfg =[];
cfg.method = 'summary';
cfg.alim = 1e-5;
cfg.keepchannel='no';
cfg.keeptrial='yes';
FT = ft_rejectvisual(cfg, FT);
[~, rejchan] = setdiff(orig_chans, FT.label);
[~, ~, rejtrial] = intersect(FT.cfg.artfctdef.summary.artifact(:,1),FT.sampleinfo(:,1));

if strcmp(format,'EEGLAB')
    EEG = eeg_interp(EEG, rejchan);
    EEG = pop_select(EEG, 'notrial', rejtrial);
elseif strcmp(format,'SPM')
    EEG = badchannels(EEG, rejchan, 1);
    EEG = badtrials(EEG, rejtrial, 1);
end
