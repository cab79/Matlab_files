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
        
        % then decide which pairing to use
        pair_idx = h.Settings.assoc.pair(conduni(c));
        
        % finally, decide actual stim according to which cue it is
        for t = 1:length(cond_ind{b}{c})
            % randomly distribute the two different possible 
            h.Seq.signal(2,cond_ind{b}{c}(t)) = h.Settings.assoc.pairing{pair_idx,h.Seq.signal(1,cond_ind{b}{c}(t))}(2);
        end
        
    end
    
    % diagnostics for debugging
    if 0
        for p = 1:2
            for st = 1:2
                num_stimtype{b}(p,st) = length(find(h.Seq.signal(st,block_ind{b})==p));
            end
        end
    end
    
end

