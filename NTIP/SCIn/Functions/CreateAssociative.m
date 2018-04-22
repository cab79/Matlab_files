function h = CreateAssociative(h)

stimtypes = [1 2];

% initialise h.Seq.signal
h.Seq.signal = nan(length(stimtypes),size(h.Seq.signal,2));

% get block indices
blockuni=unique(h.Seq.blocks);
for b = 1:length(blockuni)
    block_ind{b} = find(h.Seq.blocks==blockuni(b));
end


for b = 1:length(blockuni)
    
    % assign values to second stim based on condnum
    conduni=unique(h.Seq.condnum(block_ind{b}));
    
    for c = 1:length(conduni)
        % for each block, get index of trials for each condnum
        cond_ind{b}{c} = find(h.Seq.blocks==blockuni(b) & h.Seq.condnum==conduni(c));
        
        % assign random values for 1st cue
        randval = repmat(stimtypes,1,ceil(length(cond_ind{b}{c})/length(stimtypes)));
        randval = randval(randperm(length(randval)));
        h.Seq.signal(1,cond_ind{b}{c}) = randval(1:length(cond_ind{b}{c}));
        
        % then decide which pairing to use.
        % Depending on the probability during this block, one pairing may
        % be more likely than another
        pair_idx = h.Settings.assoc.pair(conduni(c));
        
        % for the number of unique values in each cell of
        % h.Settings.assoc.pairing (i.e. the number of stimulus outcomes
        % per pairing/cue), create a randomised sequence according to a
        % certain probability.
        % need to be balanced across different cue types.
        
        %initialise index
        cuestim_idx = nan(1,length(cond_ind{b}{c}));
        %for each cue type, get index of each type 
        cuetype = unique(h.Seq.signal(1,cond_ind{b}{c}));
        for ct = 1:length(cuetype)
            cue_idx{ct} = find(h.Seq.signal(1,cond_ind{b}{c})==cuetype(ct));
            % get number of stims per stimtype and cuetype
            stim_num{ct} = round(h.Settings.assoc.probstim * length(cue_idx{ct}));
            stim_idx{ct} = [];
            for i = 1:length(stim_num{ct})
                stim_idx{ct} = [stim_idx{ct} i*ones(1,stim_num{ct}(i))];
            end
            stim_idx{ct} = stim_idx{ct}(randperm(length(stim_idx{ct})));
            cuestim_idx(cue_idx{ct}) = stim_idx{ct};
        end
        
        % finally, decide actual stim according to which cue it is (i.e.
        % first row of signal)
        for t = 1:length(cond_ind{b}{c})
            h.Seq.signal(2,cond_ind{b}{c}(t)) = h.Settings.assoc.pairing{pair_idx,h.Seq.signal(1,cond_ind{b}{c}(t))}(cuestim_idx(t));
        end
        
    end
    
    % diagnostics for debugging
    if 1
        for p = 1:4 % for each value of signal (in either row)
            for st = 1:2
                num_stimtype{b}(p,st) = length(find(h.Seq.signal(st,block_ind{b})==p));
            end
        end
    end
    
end

