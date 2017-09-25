function h = AuditoryCreateSequence(h)

% create random list of indices of duropt rows
probs = h.Settings.duropt(:,3);
minprob = min(probs); % min prob is the divisor
mult = probs/minprob; % multiplier is the min number of repetitions of each duration option (rows of duropt)
tot = sum(mult); % total number of dur pairs
totdur = sum(sum(h.Settings.duropt(:,1:2),2) .* mult);% total duration of one set of dur pairs
num_sets = ceil(h.Settings.dur/totdur);% number of sets that can provide at least h.Settings.dur of stimulation
mod_dur = num_sets*totdur; % modified duration
% create non-randomised indices of a single set
setind = [];
for i = 1:length(mult)
    setind = [setind i*ones(1,mult(i))];
end

% create a different randomised list (block) for each repeat of the set
randind = [];
for i = 1:num_sets 

    % find sequence in which oddball trials are apart by at least nX standards
    nX = 2;

    % remove first nX standards - not to be randomised, but added to the
    % start of each set later
    setindnX = setind(nX+1:end);

    sequence_found = false;
    while ~sequence_found

        candidate = setindnX(randperm(length(setindnX)));

        w = [false candidate==1 false]; %// "close" v with zeros, and transform to logical
        starts = find(w(2:end) & ~w(1:end-1)); %// find starts of runs of non-zeros
        ends = find(~w(2:end) & w(1:end-1))-1; %// find ends of runs of non-zeros
        result = cell2mat(arrayfun(@(s,e) length(candidate(s:e)), starts, ends, 'uniformout', false)); %// build result

        % must also be no consequtive oddballs
        cand_odd = candidate>1;
        diffcand = [diff(cand_odd) 0];

        if all(result>=nX) && all(diffcand(cand_odd) ~= 0) %// check if no repeated values
            sequence_found = true;
        end
    end

    disp(['SETUP SEQUENCE: Set ' num2str(i) '/' num2str(num_sets) ' complete']);

    randind = [randind setind(1:nX) candidate];
end

t = transpose((1:mod_dur*h.Settings.fs)/h.Settings.fs);

% construct the player object: left
x = sin(2*pi*h.Settings.f0*t);
% construct the player object: right
y = sin(2*pi*(h.Settings.f0+h.Settings.df)*t);

% define durations - for duropt durations
i12 = [];
to = [];
for i = 1:length(randind)
    i12 = [
        i12;
        ones(round(h.Settings.duropt(randind(i),1)*h.Settings.fs),1); 
        zeros(round(h.Settings.duropt(randind(i),2)*h.Settings.fs),1)];
    if h.Settings.tone_dur>0
        to = [
            to;
            ones(round(h.Settings.tone_dur*h.Settings.fs),1); zeros(round((h.Settings.duropt(randind(i),1)-h.Settings.tone_dur)*h.Settings.fs),1);
            ones(round(h.Settings.tone_dur*h.Settings.fs),1); zeros(round((h.Settings.duropt(randind(i),2)-h.Settings.tone_dur)*h.Settings.fs),1);
        ];
    else
        to = ones(length(t),1);
    end
    disp(['SETUP SEQUENCE: Pair ' num2str(i) '/' num2str(length(randind)) ' appended']);
end

if length(i12)~=length(t)
    error('lengths do not match')
end

% pitch changes
if h.Settings.fpitch>0
    % alternate sin waves
    x2 = sin(2*pi*(h.Settings.f0+h.Settings.pitchdiff)*t);
    y2 = sin(2*pi*(h.Settings.f0+h.Settings.df+h.Settings.pitchdiff)*t);

    % define durations - for identical durations
    %i1 = ones((1/fpitch)*h.Settings.fs,1); 
    %i2 = zeros((1/fpitch)*h.Settings.fs,1); 
    %i12 = repmat([i1;i2],mod_dur*(h.Settings.fs/(length(i1)+length(i2))),1);

    % splice them in
    x(find(i12)) = x2(find(i12));
    y(find(i12)) = y2(find(i12));
end

% intensity changes
if h.Settings.finten>0
    % alternate sin waves
    x2 = h.Settings.intendiff*sin(2*pi*(h.Settings.f0)*t);
    y2 = h.Settings.intendiff*sin(2*pi*(h.Settings.f0+h.Settings.df)*t);

    % define durations - for identical durations
    %i1 = ones((1/finten)*h.Settings.fs,1); 
    %i2 = zeros((1/finten)*h.Settings.fs,1); 
    %i12 = repmat([i1;i2],mod_dur*(h.Settings.fs/(length(i1)+length(i2))),1);

    % splice them in
    x(find(i12)) = x2(find(i12));
    y(find(i12)) = y2(find(i12));
end
x = x.*to;
y = y.*to;

h.Seq.signal = audioplayer([x y], h.Settings.fs);
h.Seq.blocks = ones(1,1:length(h.Seq.signal));