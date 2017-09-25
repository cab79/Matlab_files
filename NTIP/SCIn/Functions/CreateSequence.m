function h = CreateSequence(h)

% create random list of indices of stimdur rows
probs = h.Settings.durprob;
minprob = min(probs); % min prob is the divisor
mult = probs/minprob; % multiplier is the min number of repetitions of each duration option (rows of duropt)
tot = sum(mult); % total number of dur pairs
totdur = sum(sum(h.Settings.stimdur,2) .* mult);% total duration of one set of dur pairs
num_sets = ceil(h.Settings.totdur/totdur);% number of sets that can provide at least h.Settings.dur of stimulation
h.totdur = num_sets*totdur; % modified duration
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

h.Seq.signal = randind; % type of signal for each trial: intensity, pitch, duration or channel
%h.Seq.pattern = ; % type of temporal pattern of the signal within each trial
%h.Seq.condnum = design(2,:);
%h.Seq.changedist = design(3,:);
h.Seq.blocks = ones(1,length(randind));

% create all trials if design is continuous
if isfield(h.Settings,'stimcontrol') && strcmp(h.Settings.design,'continuous') && h.Settings.savesinwave
    if ~isempty(h.Settings.stimcontrol)
        opt = 'create';
        h = stimtrain(h,opt);
    end
end
