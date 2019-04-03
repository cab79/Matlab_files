function [EEG,rejtrial,rejchan] = FTrejman(EEG,filterset,varargin)
%convert to FT and manually remove chan and trial condidering a particular freq

if ~isempty(varargin)
    chan = varargin{1};
end

if ~isempty(varargin) && length(varargin)>1
    format = varargin{2};
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
if ~exist('chan','var')
    chan = 1:length(orig_chans);
end
if any(chan>length(orig_chans))
    chan = 1:length(orig_chans);
end
cfg.channel = orig_chans(chan);
cfg.keepchannel='no';
cfg.keeptrial='yes';
FT = ft_rejectvisual(cfg, FT);
[~, rejchan] = setdiff(orig_chans(chan), FT.label);
[~, ~, rejtrial] = intersect(FT.cfg.artfctdef.summary.artifact(:,1),FT.sampleinfo(:,1));

if strcmp(format,'EEGLAB')
    EEG = eeg_interp(EEG, rejchan);
    if length(EEG.chanlocs)>EEG.nbchan
        EEG.chanlocs(rejchan)=[];
    end
    EEG = pop_select(EEG, 'notrial', rejtrial);
elseif strcmp(format,'SPM')
    EEG = badchannels(EEG, rejchan, 1);
    EEG = badtrials(EEG, rejtrial, 1);
end
