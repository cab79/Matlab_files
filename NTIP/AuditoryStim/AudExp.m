% Auditory experiment

handles.dur = 1;
handles.fs = 96000;
handles.f0 = 200; 
handles.df = 10;

% 1/handles.df*0.25 AND (1/df+1/pitchdiff)*0.25 needs to be a multiple of 1/f0 to prevent clicking caused by phase re-setting of the f0 sound
% settings that work:
% df=10Hz: f0=200, pitchdiff=200 
% df=25Hz: f0=200, pitchdiff=200 
% df=40Hz: f0=160, pitchdiff=160
%strategy: decide on df first, then pitchdiff (ideally a multiple of df),
%then calculate:
% 1./((1/handles.df+1/pitchdiff)*0.25 ./ [1:10])
% and choose a frequency for f0.
% seems to work if f0 and pitchdiff are both multiples of df, and f0 and
% pitchdiff are also multiples of each other?

% ERP / freq tag settings
fpitch = 2.5; % Hz - 1/fpitch must be an integer multiple of 1/handles.df 
pitchdiff = 200;
finten = 0; % Hz - 1/finten must be an integer multiple of 1/handles.df 
intendiff = 0.75; % multiple of normal intensity (value between 0 and 1)
tone_dur = 0; % set to 0 for tones to fill whole trial

% duration settings: all possible options to be randomised, and their
% probability
% left = 1st inten/pitch
% middle = 2nd inten/pitch
% right = probability
fc = max(finten, fpitch);
duropt = [
    % standard
    1/fc 1/fc 0.8
    % oddballs on left
    1/fc*0.5+1/handles.df*0.5, 1/fc, 0.025
    1/fc*1.5-1/handles.df*0.5, 1/fc, 0.025
    1/fc*0.5+1/handles.df*0.25, 1/fc, 0.025
    1/fc*1.5-1/handles.df*0.25, 1/fc, 0.025
    % oddballs on right
    1/fc, 1/fc*0.5+1/handles.df*0.5, 0.025
    1/fc, 1/fc*1.5-1/handles.df*0.5, 0.025
    1/fc, 1/fc*0.5+1/handles.df*0.25, 0.025
    1/fc, 1/fc*1.5-1/handles.df*0.25, 0.025
    ];

%% RUN

% create random list of indices of duropt rows
probs = duropt(:,3);
minprob = min(probs); % min prob is the divisor
mult = probs/minprob; % multiplier is the min number of repetitions of each duration option (rows of duropt)
tot = sum(mult); % total number of dur pairs
totdur = sum(sum(duropt(:,1:2),2) .* mult);% total duration of one set of dur pairs
num_sets = ceil(handles.dur/totdur);% number of sets that can provide at least handles.dur of stimulation
mod_dur = num_sets*totdur; % modified duration
% create non-randomised indices of a singel set
setind = [];
for i = 1:length(mult)
    setind = [setind i*ones(1,mult(i))];
end

% find sequence in which oddball trials are apart by at least nX standards
nX = 2;
sequence_found = false;
while ~sequence_found
    
    % create a different randomised list (block) for each repeat of the set
    randind = [];
    for i = 1:num_sets 
        randind = [randind setind(randperm(length(setind)))];
    end
    w = [false randind==1 false]; %// "close" v with zeros, and transform to logical
    starts = find(w(2:end) & ~w(1:end-1)); %// find starts of runs of non-zeros
    ends = find(~w(2:end) & w(1:end-1))-1; %// find ends of runs of non-zeros
    result = cell2mat(arrayfun(@(s,e) length(randind(s:e)), starts, ends, 'uniformout', false)); %// build result
    
    % must also be no consequtive oddballs
    cand_odd = randind>1;
    diffcand = [diff(cand_odd) 0];
    
    if all(result>=nX) && all(diffcand(cand_odd) ~= 0) %// check if no repeated values
        sequence_found = true;
    end
   
end

t = transpose((1:mod_dur*handles.fs)/handles.fs);

% construct the player object: left
x = sin(2*pi*handles.f0*t);
% construct the player object: right
y = sin(2*pi*(handles.f0+handles.df)*t);

% define durations - for duropt durations
i12 = [];
to = [];
for i = 1:length(randind)
    i12 = [
        i12;
        ones(round(duropt(randind(i),1)*handles.fs),1); 
        zeros(round(duropt(randind(i),2)*handles.fs),1)];
    if tone_dur>0
        to = [
            to;
            ones(round(tone_dur*handles.fs),1); zeros(round((duropt(randind(i),1)-tone_dur)*handles.fs),1);
            ones(round(tone_dur*handles.fs),1); zeros(round((duropt(randind(i),2)-tone_dur)*handles.fs),1);
        ];
    else
        to = ones(length(t),1);
    end
end

if length(i12)~=length(t)
    error('lengths do not match')
end

% pitch changes
if fpitch>0
    % alternate sin waves
    x2 = sin(2*pi*(handles.f0+pitchdiff)*t);
    y2 = sin(2*pi*(handles.f0+handles.df+pitchdiff)*t);
    
    % define durations - for identical durations
    %i1 = ones((1/fpitch)*handles.fs,1); 
    %i2 = zeros((1/fpitch)*handles.fs,1); 
    %i12 = repmat([i1;i2],mod_dur*(handles.fs/(length(i1)+length(i2))),1);
    
    % splice them in
    x(find(i12)) = x2(find(i12));
    y(find(i12)) = y2(find(i12));
end

% intensity changes
if finten>0
    % alternate sin waves
    x2 = intendiff*sin(2*pi*(handles.f0)*t);
    y2 = intendiff*sin(2*pi*(handles.f0+handles.df)*t);
    
    % define durations - for identical durations
    %i1 = ones((1/finten)*handles.fs,1); 
    %i2 = zeros((1/finten)*handles.fs,1); 
    %i12 = repmat([i1;i2],mod_dur*(handles.fs/(length(i1)+length(i2))),1);
    
    % splice them in
    x(find(i12)) = x2(find(i12));
    y(find(i12)) = y2(find(i12));
end
x = x.*to;
y = y.*to;

LRchan = audioplayer([x y], handles.fs);

choice = questdlg('Sequence created. Start sound?', ...
	'Choice', ...
	'Yes','Yes');

% play
play(LRchan)

choice = questdlg('Sound started. Stop?', ...
	'Choice', ...
	'Yes','Yes');
switch choice
    case 'Yes'
    % stop
    stop(LRchan)
end