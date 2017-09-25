function h = TSOT_Settings(h,opt)

% FILENAME OF SEQUENCE CREATION FUNCTION (without .m)
h.SeqFun = 'TSOT_CreateSequence';

switch opt
    case 'setoptions'
        
    % settings options
    h.SettingsOptions = {'1','2','3','4'};
    
    case '1'
        %% practice
        % set general options
        h = setgeneral(h);
        %record response, 1=yes, 0=no
        h.Settings.record_response=1;
        % require button press to start and re-start after a pause? 0 or 1
        h.Settings.buttonstart = 1;
        %between-train frequency (1000 divided by ISI in ms)
        h.Settings.freq = [1000./[1000]];
        %D188 output ports; rows = conditions (stimuli randomised together),
        %columns = output points within each condition
        %h.Settings.digiouts_init = [2 3; 1 3; 5 6; 4 6];
        %if opt.onehand
        %    h.Settings.digiouts_init = [2 3; 1 4];
        %    %number of hands being stimulated
        %    h.Settings.num_hands = 1;
        %else
            h.Settings.digiouts_init = [2 3; 1 4; 6 7; 5 8];
            %number of hands being stimulated
            h.Settings.num_hands = 2;
        %end
        %condition numbers to associate with each condition; rows = conditions,
        %column1 = no change conditions, column 2 = change condition.
        %if opt.onehand
        %    h.Settings.cond_num_init = [1 3; 2 4];
        %else
            h.Settings.cond_num_init = [1 3; 2 4; 5 7; 6 8];
        %end

        %on each hand, number of conditions
        h.Settings.num_changes = 2;
        %list of probabilities of a stimulus changing location on each trial
        h.Settings.change_prob_init = [0.3];
        %for each change probability (rows), list of number of trials in which the stimuli remain in same
        %location - these will be randomised within a block
        h.Settings.cp_stims_init = [1 3 6];
        %number of blocks of each change probability
        h.Settings.cond_rep_init = [1];
        %number of stimulus location changes within each condition
        h.Settings.nchanges_per_cond = [6]; 
        %number of blocks in which experiment is run until a pause
        h.Settings.nblocks = 1; 
        %distribute conditions equally among blocks
        h.Settings.distblocks = 0; 
            
    case '2'
        %% with responses
        % set general options
        h = setgeneral(h);
        %record response, 1=yes, 0=no
        h.Settings.record_response=1;
        % require button press to start and re-start after a pause? 0 or 1
        h.Settings.buttonstart = 1;
        %between-train frequency (1000 divided by ISI in ms)
        h.Settings.freq = [1000./[1000]];
        %D188 output ports; rows = conditions (stimuli randomised together),
        %columns = output points within each condition
        %h.Settings.digiouts_init = [2 3; 1 3; 5 6; 4 6];
        %if opt.onehand
        %    h.Settings.digiouts_init = [2 3; 1 4];
        %    %number of hands being stimulated
        %    h.Settings.num_hands = 1;
        %else
            h.Settings.digiouts_init = [2 3; 1 4; 6 7; 5 8];
            %number of hands being stimulated
            h.Settings.num_hands = 2;
        %end
        %condition numbers to associate with each condition; rows = conditions,
        %column1 = no change conditions, column 2 = change condition.
        %if opt.onehand
        %    h.Settings.cond_num_init = [1 3; 2 4];
        %else
            h.Settings.cond_num_init = [1 3; 2 4; 5 7; 6 8];
        %end
        %on each hand, number of conditions
        h.Settings.num_changes = 2;
        %list of probabilities of a stimulus changing location on each trial
        h.Settings.change_prob_init = [0.1 0.3 0.5];
        %for each change probability (rows), list of number of trials in which the stimuli remain in same
        %location - these will be randomised within a block
        h.Settings.cp_stims_init = [5 10 15; 1 3 6; 1 2 3];
        %number of blocks of each change probability
        h.Settings.cond_rep_init = [5 2 1];
        %number of stimulus location changes within each condition
        h.Settings.nchanges_per_cond = [30, 36, 30]; 
        %number of blocks in which experiment is run until a pause
        h.Settings.nblocks = 2; 
        %distribute conditions equally among blocks
        h.Settings.distblocks = 0; 
        
    case '3'
        %% practice
        % set general options
        h = setgeneral(h);
        %record response, 1=yes, 0=no
        h.Settings.record_response=1;
        % require button press to start and re-start after a pause? 0 or 1
        h.Settings.buttonstart = 0;
        %between-train frequency (1000 divided by ISI in ms)
        h.Settings.freq = [1000./[400]];
        %D188 output ports; rows = conditions (stimuli randomised together),
        %columns = output points within each condition
        %h.Settings.digiouts_init = [2 3; 1 3; 5 6; 4 6];
        %if opt.onehand
        %    h.Settings.digiouts_init = [2 3; 1 4];
        %    %number of hands being stimulated
        %    h.Settings.num_hands = 1;
        %else
            h.Settings.digiouts_init = [2 3; 1 4; 6 7; 5 8];
            %number of hands being stimulated
            h.Settings.num_hands = 2;
        %end
        %condition numbers to associate with each condition; rows = conditions,
        %column1 = no change conditions, column 2 = change condition.
        %if opt.onehand
        %    h.Settings.cond_num_init = [1 3; 2 4];
        %else
            h.Settings.cond_num_init = [1 3; 2 4; 5 7; 6 8];
        %end
        %on each hand, number of conditions
        h.Settings.num_changes = 2;
        %list of probabilities of a stimulus changing location on each trial
        h.Settings.change_prob_init = [0.3];
        %for each change probability (rows), list of number of trials in which the stimuli remain in same
        %location - these will be randomised within a block
        h.Settings.cp_stims_init = [1 3 6];
        %number of blocks of each change probability
        h.Settings.cond_rep_init = [1];
        %number of stimulus location changes within each condition
        h.Settings.nchanges_per_cond = [24]; 
        %number of blocks in which experiment is run until a pause
        h.Settings.nblocks = 1; 
        %distribute conditions equally among blocks
        h.Settings.distblocks = 0; 
        
    case '4'
        %% without responses
        % set general options
        h = setgeneral(h);
        %record response, 1=yes, 0=no
        h.Settings.record_response=0;
        % require button press to start and re-start after a pause? 0 or 1
        h.Settings.buttonstart = 0;
        %between-train frequency (1000 divided by ISI in ms)
        h.Settings.freq = [1000./[400]];
        %D188 output ports; rows = conditions (stimuli randomised together),
        %columns = output points within each condition
        %h.Settings.digiouts_init = [2 3; 1 3; 5 6; 4 6];
        %if opt.onehand
        %    h.Settings.digiouts_init = [2 3; 1 4];
        %    %number of hands being stimulated
        %    h.Settings.num_hands = 1;
        %else
            h.Settings.digiouts_init = [2 3; 1 4; 6 7; 5 8];
            %number of hands being stimulated
            h.Settings.num_hands = 2;
        %end
        %condition numbers to associate with each condition; rows = conditions,
        %column1 = no change conditions, column 2 = change condition.
        %if opt.onehand
        %    h.Settings.cond_num_init = [1 3; 2 4];
        %else
            h.Settings.cond_num_init = [1 3; 2 4; 5 7; 6 8];
        %end
        %on each hand, number of change conditions
        h.Settings.num_changes = 2;
        %list of probabilities of a stimulus changing location on each trial
        h.Settings.change_prob_init = [0.1 0.3 0.5];
        %for each change probability (rows), list of number of trials in which the stimuli remain in same
        %location - these will be randomised within a block
        h.Settings.cp_stims_init = [5 10 15; 1 3 6; 1 2 3];
        %number of blocks of each change probability
        h.Settings.cond_rep_init = [15 6 3];
        %number of stimulus location changes within each condition
        h.Settings.nchanges_per_cond = [180, 216, 180]; 
        %number of blocks in which experiment is run until a pause
        h.Settings.nblocks = 3; 
        %distribute conditions equally among blocks
        h.Settings.distblocks = 1; 
end

function h = setgeneral(h)

%Use labjack for controlling any equipment?
h.Settings.labjack=1;
% How to control stimulator? Options: audioplayer, labjack, spt
h.Settings.stimcontrol='labjack';
%record EEG, NS: netstation, BV: brainvision
h.Settings.record_EEG='NS';
%display RT, 1=yes, 0=no
h.Settings.displayRT=0;
% number of pulses in a train
h.Settings.npulses_train = 1;
%within-train frequency (ISI in ms)
h.Settings.p_freq=1000;
% buttonpress options: key: keyboard inputs. Blank for no button press
h.Settings.buttontype='key';
% range of keyboard presses indicating a recordable response
h.Settings.buttonopt = {'left arrow','right arrow'}; 
% adaptive staircase: meanings of the buttonopt
%h.Settings.adaptive.buttonfun = {'low','high'}; 
% adaptive staircase: corresponding signal values that would signify a
% correct answer
h.Settings.adaptive.signalval = [1 2];
% starting level of adaptive staircase
h.Settings.adaptive.startinglevel = 10; % arbitrary number that is used as a benchmark to the actual starting value?

%if h.Settings.num_changes*h.Settings.num_hands~=length(h.Settings.digiouts_init)
%    error('parameters incorrectly entered. Check: num_hands, num_changes');
%end
