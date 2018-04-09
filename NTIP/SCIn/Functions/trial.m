function h=trial(h,opt)%(A,freq,phase,samplerate,active,figureon)

% setup to: 
%   output either a single trials (h.Settings.trialdur>0) or multiple trials concatenated (h.Settings.trialdur=0)
%   vary across trials by intensity, duration, pitch or channel ('signal').
%   - each trial is create separately
%   vary within each trial by a pattern defined by intensity, duration, pitch or channel ('pattern')
%   - each pattern is created by a masking/splicing procedure

% https://uk.mathworks.com/help/matlab/import_export/record-and-play-audio.html#bsdl2eo-1
% "sound" simply plays the sound, "audioplayer" creates an object that
% allows puase and resume - maybe faster, but timing less important for
% fMRI.

dbstop if error

% if THRESHOLD
if isfield(h.Settings,'threshold')
    if ~isempty(h.Settings.threshold)
        h.seqtype.thresh=1;
        h.seqtype.adapt=0;
        h.seqtype.seq=0;
    end
% if ADAPTIVE
elseif isfield(h.Settings,'adaptive')
    if ~isempty(h.Settings.adaptive)
        h.seqtype.thresh=0;
        h.seqtype.adapt=1;
        h.seqtype.seq=0;
        for ad = 1:length(h.Settings.adaptive)
            h.atypes{ad} = h.Settings.adaptive(ad).type;
        end
    end
% Otherwise, using pre-programmed sequence to determine intensity
else
    h.seqtype.thresh=0;
    h.seqtype.adapt=0;
    h.seqtype.seq=1;
end

% if ODDBALL
h.seqtype.oddball=0;
if isfield(h.Settings,'oddballmethod')
    if ~isempty(h.Settings.oddballmethod)
        h.seqtype.oddball=1;
    end
end

% select channels
if isfield(h.Settings,'stimchan')
    h.chan = h.Settings.stimchan;
else
    h.chan = 1:h.Settings.nrchannels; % use all channels by default
end

% find trial(s) for which to create wave
if strcmp(h.Settings.design,'trials') && isfield(h,'i') % if a single trial has been defined in "trials" design
    h.trials = h.i;
elseif strcmp(h.Settings.design,'continuous')
    if h.Settings.ntrialsahead>0 && isfield(h,'i') % experiment is already running
        h.trials = h.i+(h.Settings.ntrialsahead-1);
    elseif h.Settings.ntrialsahead>0 && ~isfield(h,'i') % buffer needs pre-filling prior to experiment starting
        h.trials = 1:h.Settings.ntrialsahead;
    else
        h.trials = 1:length(h.Seq.signal); % otherwise concatenate all trials
    end
end 

% trial loop
h.Seq.stimseq = [];
% sample ending each trial
%h.Seq.trialend = [];
% instantaneous phase at end of trial
if ~isfield(h,'iphase')
    h.iphase = [];
end

% set default to phase alignment
if ~isfield(h.Settings,'alignphase')
    h.alignphase = 1;
elseif isempty(h.Settings.alignphase)
    h.alignphase = 1;
else
    h.alignphase = 0;
end

for tr = h.trials
    h.tr=tr;
    
    switch opt
        case 'set'
            h = trial_set_param(h);
        case 'set&construct'
            h = trial_set_param(h);
            h = contruct_wave(h);
    end
    
end

