function h = CreateSequence(h)

disp('Creating sequence...');
h.Seq =struct;

% create random list of indices of stimdur rows
probs = h.Settings.oddprob;
minprob = min(probs); % min prob is the divisor
mult = h.Settings.n_set*probs/minprob; % multiplier is the min number of repetitions of each option (row)
tot = sum(mult); % total number of dur pairs
if strcmp(h.Settings.oddballmethod,'duration')
    stimdur = h.Settings.oddballvalue;
else
    stimdur = h.Settings.stimdur;
end
totdur = sum(sum(stimdur,2) .* mult);% total duration of one set of dur pairs
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
    nX = h.Settings.sep_odd;

    % remove first nR standards - not to be randomised, but added to the
    % start of each set later
    nR = h.Settings.std_lead;
    setindnX = setind(nR+1:end);

    sequence_found = false;
    while ~sequence_found

        candidate = setindnX(randperm(length(setindnX)));

        w = [false candidate==h.Settings.standardind false]; %// "close" w with zeros, and transform to logical
        starts = find(w(2:end) & ~w(1:end-1)); %// find starts of runs of non-zeros
        ends = find(~w(2:end) & w(1:end-1))-1; %// find ends of runs of non-zeros
        result = cell2mat(arrayfun(@(s,e) length(candidate(s:e)), starts, ends, 'uniformout', false)); %// build result

        % must also be no consequtive oddballs
        if nX>0
            no_conseq=0;
            cand_odd = candidate>1;
            diffcand = [diff(cand_odd) 0];
            if all(diffcand(cand_odd) ~= 0) %// check if no repeated values
                no_conseq=1;
            end
        else
            no_conseq=1;
        end
        
        if all(result>=nX) && no_conseq 
            sequence_found = true;
        end
    end

    disp(['SETUP SEQUENCE: Set ' num2str(i) '/' num2str(num_sets) ' complete']);

    randind = [randind setind(1:nR) candidate];
end

% for roving oddball, each oddball stimulus is a persistent change in the
% stimulus characteristic (until the next oddball when it changes back).
% Hence, the stimulus types need updating here to be consistent with this:
if strcmp(h.Settings.oddballtype,'roving')
    for i = 1:length(randind)
        if i==1; 
            condnum(i) = 0; % first stimulus in a run should be indicated with 0
        elseif i==2 % initialise values
            if randind(i)==h.Settings.standardind % standard
                condnum(i) = h.Settings.standardind;
            elseif randind(i)==h.Settings.oddind % oddball
                condnum(i) = h.Settings.oddind;
            end
        elseif i>2 % values now depend on the direction of change
            if randind(i)==h.Settings.standardind % standard
                if condnum(i-1)==1
                    condnum(i) = 1;
                elseif condnum(i-1)==2
                    condnum(i) = 3;
                elseif condnum(i-1)==3
                    condnum(i) = 3;
                elseif condnum(i-1)==4
                    condnum(i) = 1;
                end
            elseif randind(i)==h.Settings.oddind % oddball
                if condnum(i-1)==1
                    condnum(i) = 2;
                elseif condnum(i-1)==2
                    condnum(i) = 4;
                elseif condnum(i-1)==3
                    condnum(i) = 4;
                elseif condnum(i-1)==4
                    condnum(i) = 2;
                end
            end
        end
    end
    stimtype = nan(1,length(condnum));
    stimtype(ismember(condnum,[0])) = 1;
    stimtype(ismember(condnum,[1 4])) = 1;
    stimtype(ismember(condnum,[2 3])) = 2;
else
    stimtype=randind;
    condnum=randind;
end

if strcmp(h.Settings.oddballtype,'roving')
    % identify number of standards and oddballs
    cn=unique(condnum);
    stan_ind = find(mod(cn,2)); % odd numbered
    odd_ind = find(~mod(cn,2)); % even numbered
    j=0;
    for i = cn
        j=j+1;
        ni(j) = length(find(condnum==i));
    end
    h.Seq.nodd = ni(odd_ind);
    h.Seq.nstan = ni(stan_ind);

    % store indices
    h.Seq.stan_ind = cn(stan_ind); % odd numbered
    h.Seq.odd_ind = cn(odd_ind); % even numbered

    % find a find a new sequence if nstan are too different; otherwise complete
    if max(h.Seq.nstan)-min(h.Seq.nstan) > 1
        h = CreateSequence(h);
    end
end

if ~isfield(h.Seq,'signal')
    h.Seq.signal = stimtype; % type of signal for each trial: intensity, pitch, duration or channel
    %h.Seq.pattern = ; % type of temporal pattern of the signal within each trial
    h.Seq.condnum = condnum;
    %h.Seq.changedist = design(3,:);
    h.Seq.blocks = ones(1,length(stimtype));


    % create all trials if design is continuous
    if isfield(h.Settings,'stimcontrol') && strcmp(h.Settings.design,'continuous') && h.Settings.savesinwave
        if ~isempty(h.Settings.stimcontrol)
            opt = 'create';
            h = stimtrain(h,opt);
        end
    end
end