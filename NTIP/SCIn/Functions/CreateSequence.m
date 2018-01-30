function h = CreateSequence(h)

disp('Creating sequence...');
h.Seq =struct;

%% define orthgonal conditions
conds = [];
ncond = 0;
nonprobcond = [];
if isfield(h.Settings,'conds')
    if ~isempty(h.Settings.conds)
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
    end
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
if nCP==0
    mult=1;
else
    mult = probs./repmat(gd,1,size(probs,2)); % multiplier is the min number of repetitions of each option (column) for each row (CP)
end

% adjust for minimum number of oddballs per set
minmult = ceil(h.Settings.n_odd_set./min(mult'))';
mult = mult.* repmat(minmult,1,size(mult,2));

% increase set size
%mult = repmat(h.Settings.n_set',1,size(mult,2)).*mult;

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
if isempty(h.Settings.n_set)
    try
        num_sets = h.Settings.n_odd./min(mult);
    catch
        num_sets = h.Settings.n_odd./min(mult');
    end
else
    num_sets = h.Settings.n_set;
end

h.totdur = sum(max(dursum,h.Settings.trialdur) * sum(num_sets*mult));% total duration of one set of all stim types

%% create sequence sets

% create non-randomised indices of a single set, separating each CP
% condition
setx = []; cpx=[]; setnum=0;
if nCP>0
    for cp = 1:nCP
        setind{cp} = [];
        for ii = 1:length(mult(cp,:))
            setind{cp} = [setind{cp} ii*ones(1,mult(cp,ii))];
        end

        % create a different randomised list (block) for each repeat of the set
        randind{cp} = [];
        for s = 1:num_sets(cp)

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

            disp(['SETUP SEQUENCE: CP ' num2str(cp) '/' num2str(nCP) ', Set ' num2str(s) '/' num2str(num_sets(cp)) ' complete']);

            randind{cp}{s} = [setind{cp}(1:nR) candidate];
            setnum = setnum+1;
            setx = [setx setnum*ones(1,length(candidate))];
            cpx = [cpx cp*ones(1,length(candidate))]; 
        end

        % for roving oddball, each oddball stimulus is a persistent change in the
        % stimulus characteristic (until the next oddball when it changes back).
        % Hence, the stimulus types need updating here to be consistent with this:
        
        for s = 1:num_sets(cp)
            if strcmp(h.Settings.oddballtype,'roving')
                for i = 1:length(randind{cp}{s})
                    if i==1; 
                        condnum{cp}{s}(i) = 0; % first stimulus in a run should be indicated with 0
                    elseif i==2 % initialise values
                        if randind{cp}{s}(i)==h.Settings.standardind % standard
                            condnum{cp}{s}(i) = h.Settings.standardind;
                        elseif randind{cp}{s}(i)==h.Settings.oddind % oddball
                            condnum{cp}{s}(i) = h.Settings.oddind;
                        end
                    elseif i>2 % values now depend on the direction of change
                        if randind{cp}{s}(i)==h.Settings.standardind % standard
                            if condnum{cp}{s}(i-1)==1
                                condnum{cp}{s}(i) = 1;
                            elseif condnum{cp}{s}(i-1)==2
                                condnum{cp}{s}(i) = 3;
                            elseif condnum{cp}{s}(i-1)==3
                                condnum{cp}{s}(i) = 3;
                            elseif condnum{cp}{s}(i-1)==4
                                condnum{cp}{s}(i) = 1;
                            end
                        elseif randind{cp}{s}(i)==h.Settings.oddind % oddball
                            if condnum{cp}{s}(i-1)==1
                                condnum{cp}{s}(i) = 2;
                            elseif condnum{cp}{s}(i-1)==2
                                condnum{cp}{s}(i) = 4;
                            elseif condnum{cp}{s}(i-1)==3
                                condnum{cp}{s}(i) = 4;
                            elseif condnum{cp}{s}(i-1)==4
                                condnum{cp}{s}(i) = 2;
                            end
                        end
                    end
                end
                
                stimtype{cp}{s} = nan(1,length(condnum{cp}{s}));
                stimtype{cp}{s}(ismember(condnum{cp}{s},[0])) = 1;
                stimtype{cp}{s}(ismember(condnum{cp}{s},[1 4])) = 1;
                stimtype{cp}{s}(ismember(condnum{cp}{s},[2 3])) = 2;
            else
                stan_odd_val = [1 2];
                if all(ismember(unique(randind{cp}{s}),stan_odd_val));
                    %update values according to oddball method
                    if strcmp(h.Settings.oddballmethod,'intensityindex')
                        if ~isempty(h.Settings.oddballvalue)
                            stan_odd_val = h.Settings.oddballvalue{cp};
                        end
                    end
                    stimtype{cp}{s}(randind{cp}{s}==1)=stan_odd_val(1);
                    stimtype{cp}{s}(randind{cp}{s}==2)=stan_odd_val(2);
                    condnum{cp}{s}=stimtype{cp}{s};
                else
                    stimtype{cp}{s}=randind{cp}{s};
                    condnum{cp}{s}=randind{cp}{s};
                end
            end
        end

        if strcmp(h.Settings.oddballtype,'roving')
            condnum_allset = cat(2,condnum{cp}{:});
            % identify number of standards and oddballs
            cn=unique(condnum_allset);
            stan_ind = find(mod(cn,2)); % odd numbered
            odd_ind = find(~mod(cn,2)); % even numbered
            j=0;
            for i = cn
                j=j+1;
                ni(j) = length(find(condnum_allset==i));
            end
            % finds the numbers of two types of oddballs: those occurring
            % in change from stim1 to stim2, and those vice versa.
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
            for s = 1:num_sets(cp)
                % find max value of previous CP condition
                maxval = max(condnum{cp-1}{s});
                % add this value to the condnum
                if strcmp(h.Settings.oddballtype,'classical')
                    condnum{cp}{s} = condnum{cp}{s}+maxval; 
                elseif strcmp(h.Settings.oddballtype,'roving')
                    condnum{cp}{s}(2:end) = condnum{cp}{s}(2:end)+maxval; % keep first value as 0
                end
            end
        end
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
    
    %randomise sets
    if any(h.Settings.rand_set)
        if length(h.Settings.rand_set)>1
            % get set indices to randomise
            setrand = find(ismember(cpx,find(h.Settings.rand_set)));
        else
            setrand = 1:length(setx);
        end
        % get set indices not beign randomised
        setnorand = 1:length(setx);
        setnorand(setrand)=[];
        
        % create set index
        us = unique(setx(setrand));
        us_pos = ismember(unique(setx),us);
        usnr = unique(setx(setnorand));
        usnr_pos = ismember(unique(setx),usnr);
        randus = us(randperm(length(us)));
        randnorand(us_pos) = randus;
        randnorand(usnr_pos) = usnr;
        setx_ind = [];
        for s = 1:length(randnorand)
            setx_ind = [setx_ind setx(setx==randnorand(s))];
        end
    else
        setx_ind=setx;
    end
    
    %if nCP>0
        cps=0;
        for cp = 1:nCP
            for s = 1:num_sets(cp)
                cps=cps+1;
                h.Seq.signal(setx_ind==cps) = stimtype{cp}{s}; % type of signal for each trial: intensity, pitch, duration or channel
                %h.Seq.pattern = ; % type of temporal pattern of the signal within each trial
                h.Seq.condnum(setx_ind==cps) = condnum{cp}{s};
                
                % if we are working on a cp condition whose sets are
                % randomised, balance the conds between blocks:
                randcp = ismember(cp,find(h.Settings.rand_set));
                nrandcp = length(randcp);
                if any(randcp)
                    if strcmp(h.Settings.blockopt,'cond')
                        h.Seq.blocks(setx_ind==cps) = cp*ones(1,length(stimtype{cp}{s}));
                    elseif strcmp(h.Settings.blockopt,'divide')
                        bv = s+(nCP-nrandcp); % block value
                        h.Seq.blocks(setx_ind==cps) = bv*ones(1,length(stimtype{cp}{s}));
                    end
                % otherwise assign block value according to CP value
                % not ideal - current only works if non-rand CP==1
                else
                    h.Seq.blocks(setx_ind==cps) = cp*ones(1,length(stimtype{cp}{s}));
                end
            end
        end
        [~,blockind] = sort(h.Seq.blocks);
        h.Seq.blocks = h.Seq.blocks(blockind);
        h.Seq.signal = h.Seq.signal(blockind);
        h.Seq.condnum = h.Seq.condnum(blockind);
    %else
    %    h.Seq.signal = ones(1,num_sets); % type of signal for each trial: intensity, pitch, duration or channel
    %    %h.Seq.pattern = ; % type of temporal pattern of the signal within each trial
    %    h.Seq.condnum = ones(1,num_sets);
    %    %h.Seq.changedist = design(3,:);
    %    if strcmp(h.Settings.blockopt,'cond')
    %        h.Seq.blocks = ones(1,num_sets);
    %    end
    %end
    
    % create all trials if design is continuous
    if isfield(h.Settings,'stimcontrol') && strcmp(h.Settings.design,'continuous') && h.Settings.savesinwave
        if ~isempty(h.Settings.stimcontrol)
            opt = 'create';
            h = stimtrain(h,opt);
        end
    end
    figure;plot(h.Seq.condnum)
end

%% create Adaptive type order

if isfield(h.Settings,'adaptive')
    if length(h.Settings.adaptive)>1 && isfield(h.Settings,'adaptive_general')
        
        % Create order of runs
        if strcmp(h.Settings.adaptive_general.seqtype,'alt')
            order = repmat(h.Settings.adaptive_general.adapttypes,1,h.Settings.adaptive(1).nRuns);
        elseif strcmp(h.Settings.adaptive_general.seqtype,'rand')
            order=[];
            nRuns=[h.Settings.adaptive.nRuns];
            bs=h.Settings.adaptive_general.seqrandblocksize;
            num_blocks = sum(nRuns)/bs;
            for nb=1:num_blocks
                bi = bs*(nb-1)+1:bs*nb;
                runs = round(([nRuns(1)/sum(nRuns),nRuns(2)/sum(nRuns)] * bs));
                order_bi=[];
                for i = 1:length(h.Settings.adaptive_general.adapttypes)
                    order_bi = [order_bi h.Settings.adaptive_general.adapttypes(i)*ones(1,runs(i))];
                end
                order(bi) = order_bi(randperm(length(order_bi)));
            end
        end
        
        % Create trial order
        uorder = unique(order)';
        for i = 1:length(uorder)
            ntype(i) = h.Settings.adaptive(i).trialsperrun;
        end
        
        % create sequence
        aseq=[];
        for i = 1:length(order)
            aseq = [aseq order(i)*ones(1,ntype(order(i)))];
        end
        
        % get condition numbers and their indices
        condval=[];
        if isfield(h.Settings.adaptive_general.selectcond,'cp') % use cp condition
            for cpi = 1:length(h.Settings.adaptive_general.selectcond.cp)
                condval = [condval unique(condnum{h.Settings.adaptive_general.selectcond.cp(cpi)}{1})];
            end
        end
        if ~isempty(condval)
            condind = find(ismember(h.Seq.condnum,condval));
        else
            condind = 1:length(h.Seq.condnum);
        end
        
        % check Seq.signal is long enough
        if length(h.Seq.signal)<length(condind)
            error('sequence is not long enough for this number of adaptive trials')
        end
        if length(aseq)<length(condind)
            error('number of adaptive trials is not long enough for this sequence')
        end
        
        % add to Seq
        h.Seq.adapttype = nan(1,length(h.Seq.condnum));
        h.Seq.adapttype(condind) = aseq(1:length(condind));
        
    end
    hold on; plot(h.Seq.adapttype,'r')
end
