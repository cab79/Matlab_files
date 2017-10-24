function h = OddballCreateSequence(h)

disp('Creating sequence...');
h.Seq =struct;

% stimtype = e.g. stimtypets (if multiple stimtypets stimulated)

if h.Settings.distblocks==1 && h.Settings.nblocks>1
    rep_design = h.Settings.nblocks; % create uniquely randomsised sequences for each block that are balanced by conditions
else
    rep_design = 1;
end

combdesign=[];

for rp = 1:rep_design

    %% initiate parameters (from parameters script)
    nfreq = length(h.Settings.freq); % number of frequencies (ISIs) of stimulation to use
    num_prob = size(h.Settings.change_prob_init,2); % number of probabilities of stimulus change to use
    change_prob_init = repmat(h.Settings.change_prob_init,1,h.Settings.num_oddballtypes*h.Settings.num_chantypes); % change probs for each block type
    cond_rep_init = repmat(h.Settings.cond_rep_init/rep_design,1,h.Settings.num_oddballtypes*h.Settings.num_chantypes); % block repetitions for each block type
    stimtypeouts_init = reshape(reshape(repmat(h.Settings.stimtypeouts_init,1,num_prob)',h.Settings.num_oddballtypes*h.Settings.num_chantypes*num_prob*2,1)',2,h.Settings.num_oddballtypes*h.Settings.num_chantypes*num_prob);
    cond_num_init = reshape(reshape(repmat(h.Settings.cond_num_init,1,num_prob)',h.Settings.num_oddballtypes*h.Settings.num_chantypes*num_prob*2,1)',2,h.Settings.num_oddballtypes*h.Settings.num_chantypes*num_prob);
    cp_stims_init = repmat(h.Settings.cp_stims_init',1,h.Settings.num_chantypes*h.Settings.num_oddballtypes); % all possible N stimuli occuring before an oddball change

    % if there are a number of probabilities (of stimulus change) during the
    % experiment, we need separate condition numbers for each.
    if num_prob>1
        add_cond_num = max(cond_num_init(:));
        for c = 1:size(cond_num_init,2)
            cind = find(cond_num_init(1,c)==unique(cond_num_init(1,:)));
            d=c-(cind-1)*num_prob;
            cond_num_init(:,c) = cond_num_init(:,c) + (d-1)*add_cond_num*ones(size(cond_num_init(:,c)));
        end
    end

    cond_rep = [];
    change_prob = [];
    stimtypeouts = [];
    cond_num = [];
    cp_stims = [];
    for i = 1:length(cond_rep_init) % for each condition
        change_prob_temp = [];
        cond_rep_temp = [];
        stimtypeouts_temp = [];
        cond_num_temp = [];
        cp_stims_temp = [];
        for j = 1:cond_rep_init(i) % for each repetition of that condition
            cond_rep_temp = [cond_rep_temp cond_rep_init(:,i)];
            change_prob_temp = [change_prob_temp change_prob_init(:,i)];
            stimtypeouts_temp = [stimtypeouts_temp stimtypeouts_init(:,i)];
            cp_stims_temp = [cp_stims_temp cp_stims_init(:,i)];
            cond_num_temp = [cond_num_temp cond_num_init(:,i)];
        end
        cond_rep = [cond_rep cond_rep_temp];
        change_prob = [change_prob change_prob_temp];
        stimtypeouts = [stimtypeouts stimtypeouts_temp];
        cp_stims = [cp_stims cp_stims_temp];
        cond_num = [cond_num cond_num_temp];
    end

    [~,randseq] = sort(rand([1,sum(cond_rep_init)]));

    change_prob = change_prob(:,randseq);
    cond_rep = cond_rep(:,randseq);
    stimtypeouts = stimtypeouts(:,randseq)';
    cp_stims = cp_stims(:,randseq);
    cond_num = cond_num(:,randseq)';

    num_randblocks = 1;

    %nstim_range = 2:6; % 1 stimulus that has changed location, plus mean of 3 stimuli have do not change, meaning 3:1 ratio of no-change:change
    %nstim_range = 1:2; % equal prob of staying or changing, i.e. 50% stay, 50% change (split over num of changes for that stimtype)
    %nstim_range = 1; % always changes with equal prob


    designs = [];
    % LOOP TO CREATE ODDBALL SEQUENCES FOR EACH CONDITION (AFTER CONDS HAVE BEEN RANDOMISED)
    for b = 1:size(stimtypeouts,1)

        stimtypeout = stimtypeouts(b,:); % e.g. each channel
        nstimtypes = length(stimtypeout);
        stimtype_changes = cell(1,nstimtypes);
        stimtype_start=stimtypeout(1); % stimulus type to start the sequence

        i = find(cond_rep_init==cond_rep(b)); % index of cond_rep_init for this condition
        nchanges_per_cond_b = h.Settings.nchanges_per_cond(i(1))/rep_design;

        nchanges_per_block = nchanges_per_cond_b/(cond_rep(b)*h.Settings.num_oddballtypes);

        if ~any(cell2mat(cellfun(@(x) ~isempty(x),stimtype_changes,'UniformOutput', false)))
            for nd = 1:nstimtypes
                stimtype_changes{1,nd} = stimtypeout;
                stimtype_changes{1,nd}(find(stimtype_changes{1,nd}==stimtypeout(nd)))=[];
            end
        end

        stimtype_change_diffs = cell(1,nstimtypes);
        for nd = 1:nstimtypes
            stimtype_change_diffs{1,nd} = abs(stimtypeout(nd)-stimtype_changes{1,nd});
        end
        changes = unique(cell2mat(stimtype_change_diffs));
        nchanges = length(changes);

        fchanges = [];
        denom=1;
        for nc = 1:nchanges
            fchanges(nc) = sum(cell2mat(cellfun(@(x) length(find(x==changes(nc))),stimtype_change_diffs,'UniformOutput', false)));
            denom=lcm(fchanges(nc),denom);
        end
        nstim_per_change = ceil(nchanges_per_block/denom)*denom; % number of stimuli for each location change (1 to 4) must be integer-divisible by the frequency of changes, so that an integer number of repeats of each change can be included

        nstim_fchange = cell(1,length(stimtype_change_diffs));
        for nd = 1:length(nstim_fchange)
            for nc = 1:length(stimtype_change_diffs{1,nd})
                stimtype_change_diffs{1,nd}(1,nc);
                nstim_fchange{1,nd}(1,nc) = nstim_per_change/fchanges(find(changes==stimtype_change_diffs{1,nd}(1,nc))); % generates the integer number of changes of each type (1 to 4)
            end
        end

        nstim_fchange_orig = nstim_fchange;

        stimtype_lists = cell(1,nstimtypes);
        for nd = 1:nstimtypes
            stimtype_lists{1,nd}=[];
            for nb = 1:num_randblocks
                block_list=[];
                for ch = 1:length(stimtype_changes{1,nd})
                    numrep = min(nstim_fchange{1,nd}(1,ch), round(nstim_fchange_orig{1,nd}(1,ch)/num_randblocks));
                    nstim_fchange{1,nd}(1,ch) = nstim_fchange{1,nd}(1,ch)-numrep;
                    block_list= [block_list repmat(stimtype_changes{1,nd}(1,ch),1,numrep)];
                end
                block_list=block_list(randperm(length(block_list)));
                stimtype_lists{1,nd} = [stimtype_lists{1,nd} block_list];
            end
        end

        nstim_range = cp_stims(:,b)';
        nstim_range(nstim_range==0)=[];

        num_stim = cell(1,nchanges);
        for nc = 1:nchanges
            stimlist = [];
            for nt = 1:nstim_per_change/length(nstim_range) % this may result in fewer stims than in stimtype_lists if nstim_per_change is not integer divisible by the range of nstim
                if length(nstim_range)==1 
                    stimlist = [stimlist nstim_range];
                else
                    stimlist = [stimlist randsample(nstim_range,length(nstim_range))];
                end
            end
            num_stim{1,nc} = stimlist;
        end

        design = [];
        for i = 1:length(num_stim{1,1})*nchanges
            if i==1
                nstim = ceil(median(nstim_range-1));
                stimtype=stimtype_start;
            elseif mean(cell2mat(cellfun(@(x) length(x),num_stim,'UniformOutput', false))) < 2 && mean(cell2mat(cellfun(@(x) length(x),num_stim,'UniformOutput', false))) > 2*std(cell2mat(cellfun(@(x) length(x),num_stim,'UniformOutput', false)))
                %next_iter==0;
                break
            else
                % num_empty = length(find(cell2mat(cellfun(@(x) isempty(x),num_stim,'UniformOutput', false))));
                if isempty(num_stim{1,find(changes==abs(stimtype-next_stimtype))}) || isempty(stimtype_lists{1,find(stimtypeout==next_stimtype)})
                    next_iter=1;
                    iteration = iteration+1
                    break
                end
                nstim = num_stim{1,find(changes==abs(stimtype-next_stimtype))}(1);
                num_stim{1,find(changes==abs(stimtype-next_stimtype))}(1) = [];
                stimtype=next_stimtype;
            end
            design = [design repmat(stimtype,1,nstim)];

            next_stimtype=stimtype_lists{1,find(stimtypeout==stimtype)}(1);
            stimtype_lists{1,find(stimtypeout==stimtype)}(1)=[];
        end

        % generate conditions
        for ns = 1:size(design,2);
            if ns==1
                design(2,ns) = 0;
                design(3,ns) = 0;
            elseif (design(1,ns)-design(1,ns-1))~=0
                design(2,ns) = cond_num(b,1); 
                design(3,ns) = abs(design(1,ns)-design(1,ns-1));
            elseif design(2,ns-1)==0
                design(2,ns) = cond_num(b,2);
                design(3,ns) = abs(design(1,ns)-design(1,ns-1));
            elseif design(2,ns-1)==cond_num(b,1);
                design(2,ns) = cond_num(b,2); 
                design(3,ns) = abs(design(1,ns)-design(1,ns-1));
            else
                design(2,ns) = design(2,ns-1);
                design(3,ns) = abs(design(1,ns)-design(1,ns-1));
            end
            if rep_design>1; design(4,ns) = rp; end;
        end

        designs = [designs design];

    end
    combdesign = [combdesign designs];
end
design = combdesign;

% if blocking does not require even distribution of conditions in each
% block:
if strcmp(h.Settings.blockopt,'divide') && rep_design==1 && h.Settings.nblocks>1
    % for each block
    for nb = 1:h.Settings.nblocks
        % find the MINIMUM trial index in the nb block
        b_end = floor(nb*(size(design,2)/h.Settings.nblocks));
        % start counting from the minimum to find the end of the nb block
        ns=b_end;

        % if not the last trial, start counting up...
        while ns < size(design,2)
            % if the first trial of a new condition, this is a good place
            % for a new block
            if design(2,ns)==0
                % 
                bn(nb)=ns; % bn = index of last trial of block nb
                break;
            end
            ns=ns+1;
        end
        % if the last trial...
        if ns == size(design,2)
            bn(nb)=ns;
        end

        % assign trial indices of each block
        if nb==1; 
            ind = 1:bn(nb);
        else
            ind = bn(nb-1)+1:bn(nb);
        end
        design(4,ind) = nb*ones(length(ind),1);
    end
elseif strcmp(h.Settings.blockopt,'cond')
    nb=0;
    for ns = 1:size(design,2)
        if design(2,ns)==0
            nb=nb+1;
        end
        design(4,ns) = nb;
    end
else
    design(4,:) = 1;
end

h.Seq.signal = design(1,:);
h.Seq.condnum = design(2,:);
h.Seq.changedist = design(3,:);
h.Seq.blocks = design(4,:);

% response probe sequence
%augment conditions
if isfield(h.Settings,'RPconds')
    aug_factor = (numel(cond_num_init)/numel(h.Settings.cond_num_init)-1) * max(h.Settings.cond_num_init(:));
    RPconds = [h.Settings.RPconds h.Settings.RPconds+aug_factor];
    RPi = find(ismember(h.Seq.condnum,RPconds)); % indices of all possible stims to use
    nSelect = round(h.Settings.RPprob*length(RPi));
    RPrand = RPi(randperm(length(RPi)));
    RPselect = RPrand(1:nSelect);
    h.Seq.RP = zeros(1,length(h.Seq.condnum));
    h.Seq.RP(RPselect) = 1;
end
% create all trials if design is continuous
if isfield(h.Settings,'stimcontrol') && strcmp(h.Settings.design,'continuous') && h.Settings.savesinwave
    if ~isempty(h.Settings.stimcontrol)
        opt = 'create';
        h = stimtrain(h,opt);
    end
end