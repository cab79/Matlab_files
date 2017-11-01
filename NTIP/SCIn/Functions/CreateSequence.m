function h = CreateSequence(h)

disp('Creating sequence...');
h.Seq =struct;

%% define orthgonal conditions
if isfield(h.Settings,'conds')
    conds = h.Settings.conds;
    allcondind='';
    for i = 1:length(conds)
        condval{i} = h.Settings.(conds{i});
        if size(condval{i},1)<2
            continue
        end
        condind{i} = 1:size(condval{i},1);
        allcondind = [allcondind ' condind{' num2str(i) '},'];
    end
    eval(['allcond = allcomb(' allcondind(1:end-1) ');']);
    allcond = 1:i;
    probcond = find(strcmp(conds,h.Settings.oddprob));
    nonprobcond = allcond; nonprobcond(probcond) = [];
else
    conds = [];
    ncond = 0;
    nonprobcond = [];
end

%% work out parameters of each sequence set (length, number of oddballs, probability of oddballs, etc.)

% ensure probs are integers for gcd calculation
probs = h.Settings.oddprob*1000000;

% number of CP conditions
nCP = size(probs,1);

% find greatest common divisor to find least number of repetitions of each
% oddball
for i = 1:nCP
    gd(i,1) = probs(i,1);
    for ii=2:numel(probs(i,:))
        gd(i,1) = gcd(gd(i,1),probs(i,ii));
    end
end
mult = probs./repmat(gd,1,size(probs,2)); % multiplier is the min number of repetitions of each option (row) for each column (CP)

% adjust for minimum number of oddballs per set
minmult = ceil(h.Settings.n_odd_set./min(mult'))';
mult = mult.* repmat(minmult,1,size(mult,2));

% increase set size
mult = repmat(h.Settings.n_set',1,size(mult,2)).*mult;

% calculate total duration of one set
%tot = sum(mult(:)); % total number of the minimum set of stim types
if strcmp(h.Settings.oddballmethod,'duration')
    stimdur = h.Settings.oddballvalue;
else
    stimdur = h.Settings.stimdur;
end
if iscell(stimdur)
    dursum = cellfun(@sum,stimdur);
else
    dursum = sum(stimdur,2);
end
totdur = sum(max(dursum,h.Settings.trialdur) .* mult(:));% total duration of one set of all stim types

% calculate number of sets that can provide at least h.Settings.totdur of stimulation AND n_odd oddballs per CP
num_sets=0;
if isfield(h.Settings,'totdur')
    num_sets = ceil(h.Settings.totdur/totdur);
end
num_sets = max(num_sets,max(ceil(h.Settings.n_odd./min(mult))));

h.totdur = num_sets*totdur; % modified duration

%% create sequence sets

% create non-randomised indices of a single set, separating each CP
% condition
for cp = 1:nCP
    setind{cp} = [];
    for ii = 1:length(mult(cp,:))
        setind{cp} = [setind{cp} ii*ones(1,mult(cp,ii))];
    end

    % create a different randomised list (block) for each repeat of the set
    randind{cp} = [];
    for i = 1:num_sets

        % find sequence in which oddball trials are apart by at least nX standards
        nX = h.Settings.sep_odd(cp);
        nXi = h.Settings.sep_odd_ind;

        % remove first nR standards - not to be randomised, but added to the
        % start of each set later
        nR = h.Settings.std_lead(cp);
        setindnX = setind{cp}(nR+1:end);

        sequence_found = false;
        while ~sequence_found

            candidate = setindnX(randperm(length(setindnX)));
            
            no_conseq=1;min_stan=1;

            for ii = 1:length(nXi)
                % find indices of oddball types NOT to consider and
                % make them 1 (as if standards)
                sub_cand = candidate;
                oi = nXi{ii};
                oi_ind = ~ismember(candidate,oi);
                sub_cand(oi_ind) = 1;
                    
                % How long is each sequence of standards?
                w = [false sub_cand==h.Settings.standardind false]; %// "close" w with zeros, and transform to logical
                starts = find(w(2:end) & ~w(1:end-1)); %// find starts of runs of non-zeros
                ends = find(~w(2:end) & w(1:end-1))-1; %// find ends of runs of non-zeros
                result = cell2mat(arrayfun(@(s,e) length(sub_cand(s:e)), starts, ends, 'uniformout', false)); %// build result
                if ~all(result>=nX)
                    min_stan=0;
                end

                % must also be no consecutive oddballs
                if nX>0
                    cand_odd = sub_cand>1;
                    diffcand = [diff(cand_odd) 0];
                    if ~all(diffcand(cand_odd) ~= 0) %// check if no repeated values
                        no_conseq=0;
                    end
                end
            end

            if min_stan && no_conseq 
                sequence_found = true;
            end
        end

        disp(['SETUP SEQUENCE: Set ' num2str(i) '/' num2str(num_sets) ' complete']);

        randind{cp} = [randind{cp} setind{cp}(1:nR) candidate];
    end

    % for roving oddball, each oddball stimulus is a persistent change in the
    % stimulus characteristic (until the next oddball when it changes back).
    % Hence, the stimulus types need updating here to be consistent with this:
    if strcmp(h.Settings.oddballtype,'roving')
        for i = 1:length(randind{cp})
            if i==1; 
                condnum{cp}(i) = 0; % first stimulus in a run should be indicated with 0
            elseif i==2 % initialise values
                if randind{cp}(i)==h.Settings.standardind % standard
                    condnum{cp}(i) = h.Settings.standardind;
                elseif randind{cp}(i)==h.Settings.oddind % oddball
                    condnum{cp}(i) = h.Settings.oddind;
                end
            elseif i>2 % values now depend on the direction of change
                if randind{cp}(i)==h.Settings.standardind % standard
                    if condnum{cp}(i-1)==1
                        condnum{cp}(i) = 1;
                    elseif condnum{cp}(i-1)==2
                        condnum{cp}(i) = 3;
                    elseif condnum{cp}(i-1)==3
                        condnum{cp}(i) = 3;
                    elseif condnum{cp}(i-1)==4
                        condnum{cp}(i) = 1;
                    end
                elseif randind{cp}(i)==h.Settings.oddind % oddball
                    if condnum{cp}(i-1)==1
                        condnum{cp}(i) = 2;
                    elseif condnum{cp}(i-1)==2
                        condnum{cp}(i) = 4;
                    elseif condnum{cp}(i-1)==3
                        condnum{cp}(i) = 4;
                    elseif condnum{cp}(i-1)==4
                        condnum{cp}(i) = 2;
                    end
                end
            end
        end
        stimtype{cp} = nan(1,length(condnum{cp}));
        stimtype{cp}(ismember(condnum{cp},[0])) = 1;
        stimtype{cp}(ismember(condnum{cp},[1 4])) = 1;
        stimtype{cp}(ismember(condnum{cp},[2 3])) = 2;
    else
        stimtype{cp}=randind{cp};
        condnum{cp}=randind{cp};
    end

    if strcmp(h.Settings.oddballtype,'roving')
        % identify number of standards and oddballs
        cn=unique(condnum{cp});
        stan_ind = find(mod(cn,2)); % odd numbered
        odd_ind = find(~mod(cn,2)); % even numbered
        j=0;
        for i = cn
            j=j+1;
            ni(j) = length(find(condnum{cp}==i));
        end
        h.Seq.nodd{cp} = ni(odd_ind);
        h.Seq.nstan{cp} = ni(stan_ind);

        % store indices
        h.Seq.stan_ind{cp} = cn(stan_ind); % odd numbered
        h.Seq.odd_ind{cp} = cn(odd_ind); % even numbered

        % find a new sequence if nstan are too different; otherwise complete
        %if max(h.Seq.nstan{cp})-min(h.Seq.nstan{cp}) > 1
        %    h = CreateSequence(h);
        %end
    end
     
    % make condition numbers distinct for different CP conditions
    if cp>1
        % find max value of previous CP condition
        maxval = max(condnum{cp-1});
        % add this value to the condnum
        condnum{cp}(2:end) = condnum{cp}(2:end)+maxval; % keep first value as 0
    end
end

%% duplicate sequences for non-probability conditions
if ~isempty(nonprobcond)
    
end

%% create final sequences/blocks
if ~isfield(h.Seq,'signal')
    h.Seq.signal=[];
    h.Seq.condnum=[];
    h.Seq.blocks=[];
    for cp = 1:nCP
        h.Seq.signal = [h.Seq.signal stimtype{cp}]; % type of signal for each trial: intensity, pitch, duration or channel
        %h.Seq.pattern = ; % type of temporal pattern of the signal within each trial
        h.Seq.condnum = [h.Seq.condnum condnum{cp}];
        %h.Seq.changedist = design(3,:);
        if strcmp(h.Settings.blockopt,'cond')
            h.Seq.blocks = [h.Seq.blocks cp*ones(1,length(stimtype{cp}))];
        end
    end
    
    % create all trials if design is continuous
    if isfield(h.Settings,'stimcontrol') && strcmp(h.Settings.design,'continuous') && h.Settings.savesinwave
        if ~isempty(h.Settings.stimcontrol)
            opt = 'create';
            h = stimtrain(h,opt);
        end
    end
end