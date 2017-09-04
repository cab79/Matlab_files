function h = TSOT_Settings(h,opt)

% FILENAME OF SEQUENCE CREATION FUNCTION (without .m)
h.SeqFun = 'TSOT_CreateSequence';

switch opt
    case 'setoptions'
        
    % settings options
    h.SettingsOptions = {'1','2','3','4'};
    
    case '1'
        %% practice
        %record response, 1=yes, 0=no
        h.record_response=1;
        %record EEG, 1=yes, 0=no
        h.record_EEG=0;
        %display RT, 1=yes, 0=no
        h.displayRT=1;
        %between-train frequency (1000 divided by ISI in ms)
        h.freq = [1000./[1000]];
        % number of pulses in a train
        h.npulses_train = 1;
        %within-train frequency (ISI in ms)
        h.p_freq=1000;
        %D188 output ports; rows = conditions (stimuli randomised together),
        %columns = output points within each condition
        %h.digiouts_init = [2 3; 1 3; 5 6; 4 6];
        %if opt.onehand
        %    h.digiouts_init = [2 3; 1 4];
        %    %number of hands being stimulated
        %    h.num_hands = 1;
        %else
            h.digiouts_init = [2 3; 1 4; 6 7; 5 8];
            %number of hands being stimulated
            h.num_hands = 2;
        %end
        %condition numbers to associate with each condition; rows = conditions,
        %column1 = no change conditions, column 2 = change condition.
        %if opt.onehand
        %    h.cond_num_init = [1 3; 2 4];
        %else
            h.cond_num_init = [1 3; 2 4; 5 7; 6 8];
        %end

        %on each hand, number of conditions
        h.num_changes = 2;
        %list of probabilities of a stimulus changing location on each trial
        h.change_prob_init = [0.3];
        %for each change probability (rows), list of number of trials in which the stimuli remain in same
        %location - these will be randomised within a block
        h.cp_stims_init = [1 3 6];
        %number of blocks of each change probability
        h.cond_rep_init = [1];
        %number of stimulus location changes within each condition
        h.nchanges_per_cond = [6]; 
        %number of blocks in which experiment is run until a pause
        h.nblocks = 1; 
        %distribute conditions equally among blocks
        h.distblocks = 0; 
            
    case '2'
        %% with responses
        %record response, 1=yes, 0=no
        h.record_response=1;
        %record EEG, 1=yes, 0=no
        h.record_EEG=1;
        %display RT, 1=yes, 0=no
        h.displayRT=1;
        %between-train frequency (1000 divided by ISI in ms)
        h.freq = [1000./[1000]];
        % number of pulses in a train
        h.npulses_train = 1;
        %within-train frequency (ISI in ms)
        h.p_freq=1000;
        %D188 output ports; rows = conditions (stimuli randomised together),
        %columns = output points within each condition
        %h.digiouts_init = [2 3; 1 3; 5 6; 4 6];
        %if opt.onehand
        %    h.digiouts_init = [2 3; 1 4];
        %    %number of hands being stimulated
        %    h.num_hands = 1;
        %else
            h.digiouts_init = [2 3; 1 4; 6 7; 5 8];
            %number of hands being stimulated
            h.num_hands = 2;
        %end
        %condition numbers to associate with each condition; rows = conditions,
        %column1 = no change conditions, column 2 = change condition.
        %if opt.onehand
        %    h.cond_num_init = [1 3; 2 4];
        %else
            h.cond_num_init = [1 3; 2 4; 5 7; 6 8];
        %end
        %on each hand, number of conditions
        h.num_changes = 2;
        %list of probabilities of a stimulus changing location on each trial
        h.change_prob_init = [0.1 0.3 0.5];
        %for each change probability (rows), list of number of trials in which the stimuli remain in same
        %location - these will be randomised within a block
        h.cp_stims_init = [5 10 15; 1 3 6; 1 2 3];
        %number of blocks of each change probability
        h.cond_rep_init = [5 2 1];
        %number of stimulus location changes within each condition
        h.nchanges_per_cond = [30, 36, 30]; 
        %number of blocks in which experiment is run until a pause
        h.nblocks = 2; 
        %distribute conditions equally among blocks
        h.distblocks = 0; 
        
    case '3'
        %% practice
        %record response, 1=yes, 0=no
        h.record_response=1;
        %record EEG, 1=yes, 0=no
        h.record_EEG=0;
        %display RT, 1=yes, 0=no
        h.displayRT=1;
        %between-train frequency (1000 divided by ISI in ms)
        h.freq = [1000./[400]];
        % number of pulses in a train
        h.npulses_train = 1;
        %within-train frequency (ISI in ms)
        h.p_freq=1000;
        %D188 output ports; rows = conditions (stimuli randomised together),
        %columns = output points within each condition
        %h.digiouts_init = [2 3; 1 3; 5 6; 4 6];
        %if opt.onehand
        %    h.digiouts_init = [2 3; 1 4];
        %    %number of hands being stimulated
        %    h.num_hands = 1;
        %else
            h.digiouts_init = [2 3; 1 4; 6 7; 5 8];
            %number of hands being stimulated
            h.num_hands = 2;
        %end
        %condition numbers to associate with each condition; rows = conditions,
        %column1 = no change conditions, column 2 = change condition.
        %if opt.onehand
        %    h.cond_num_init = [1 3; 2 4];
        %else
            h.cond_num_init = [1 3; 2 4; 5 7; 6 8];
        %end
        %on each hand, number of conditions
        h.num_changes = 2;
        %list of probabilities of a stimulus changing location on each trial
        h.change_prob_init = [0.3];
        %for each change probability (rows), list of number of trials in which the stimuli remain in same
        %location - these will be randomised within a block
        h.cp_stims_init = [1 3 6];
        %number of blocks of each change probability
        h.cond_rep_init = [1];
        %number of stimulus location changes within each condition
        h.nchanges_per_cond = [24]; 
        %number of blocks in which experiment is run until a pause
        h.nblocks = 1; 
        %distribute conditions equally among blocks
        h.distblocks = 0; 
        
    case '4'
        %% without responses
        %record response, 1=yes, 0=no
        h.record_response=0;
        %record EEG, 1=yes, 0=no
        h.record_EEG=1;
        %display RT, 1=yes, 0=no
        h.displayRT=0;
        %between-train frequency (1000 divided by ISI in ms)
        h.freq = [1000./[400]];
        % number of pulses in a train
        h.npulses_train = 1;
        %within-train frequency (ISI in ms)
        h.p_freq=1000;
        %D188 output ports; rows = conditions (stimuli randomised together),
        %columns = output points within each condition
        %h.digiouts_init = [2 3; 1 3; 5 6; 4 6];
        %if opt.onehand
        %    h.digiouts_init = [2 3; 1 4];
        %    %number of hands being stimulated
        %    h.num_hands = 1;
        %else
            h.digiouts_init = [2 3; 1 4; 6 7; 5 8];
            %number of hands being stimulated
            h.num_hands = 2;
        %end
        %condition numbers to associate with each condition; rows = conditions,
        %column1 = no change conditions, column 2 = change condition.
        %if opt.onehand
        %    h.cond_num_init = [1 3; 2 4];
        %else
            h.cond_num_init = [1 3; 2 4; 5 7; 6 8];
        %end
        %on each hand, number of change conditions
        h.num_changes = 2;
        %list of probabilities of a stimulus changing location on each trial
        h.change_prob_init = [0.1 0.3 0.5];
        %for each change probability (rows), list of number of trials in which the stimuli remain in same
        %location - these will be randomised within a block
        h.cp_stims_init = [5 10 15; 1 3 6; 1 2 3];
        %number of blocks of each change probability
        h.cond_rep_init = [15 6 3];
        %number of stimulus location changes within each condition
        h.nchanges_per_cond = [180, 216, 180]; 
        %number of blocks in which experiment is run until a pause
        h.nblocks = 3; 
        %distribute conditions equally among blocks
        h.distblocks = 1; 
end

%if h.num_changes*h.num_hands~=length(h.digiouts_init)
%    error('parameters incorrectly entered. Check: num_hands, num_changes');
%end
