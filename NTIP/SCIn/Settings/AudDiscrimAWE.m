function h = AudDiscrimAWE(h,opt)


% FILENAME OF SEQUENCE CREATION FUNCTION (without .m)
h.SeqFun = 'CreateSequence';

switch opt
    
    case 'setoptions'
        
    % settings options
    h.SettingsOptions = {'10Hz'};
    
    case '10Hz'
        
    %% TRIALS or CONTINUOUS?
    h.Settings.design = 'continuous';
    % if continuous, how many trials ahead should be in the player schedule?
    % (applied to stimulation via soundcard only)
    h.Settings.ntrialsahead = 2;  %0 = all trials
    
    %% EXPERIMENTAL CONDIITIONS
    % name the settings that define orthogonal condtions at a different row
    h.Settings.conds = {};
    
    %% Output options
    % save sinwave from all trials as part of stim sequence file
    h.Settings.savesinwave = 0;
    
    %% EQUIPMENT CONTROL
    % record EEG, NS: netstation, BV: brainvision, 'serial': serial port
    %h.Settings.record_EEG='serial';
    %h.Settings.EEGport = 'COM3'; % only needed for 'serial' EEG triggers
    %h.Settings.EEGMarkPattern = 1; % mark EEG for every change in stimulus pattern (0 = start of trial only)
    h.Settings.labjack=0; % Use labjack for controlling any equipment?
    h.Settings.stimcontrol='PsychPortAudio'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    h.Settings.nrchannels = 2; % total number of channels, e.g. on sound card
    h.Settings.stimchan = [1 2]; % channels on stimulator to use
    
    %% BLOCKING/RUN OPTIONS
    % 'divide' = equally divide trials by nblocks; 
    % 'cond' = separate block for each condition
    h.Settings.blockopt = 'cond';
    % further options for 'divide':
        % number of blocks (containing multiple conditions)
    %    h.Settings.nblocks = 2; % must integer-divide each value in h.Settings.cond_rep_init
        %distribute conditions equally among blocks
    %    h.Settings.distblocks = 1;
    % options to start sequence at beginning of every run
    % 'msgbox', 'labjack', 'buttonpress', 'audio' - can have more than one in
    % cell array
    h.Settings.blockstart = {'buttonpress'}; % audio,labjack,audio
    % names of any audiofiles
    h.Settings.audiofile = {'instruct.wav'}; % labjack
    
    %% Condition-independent stimulus parameters - can be superceded by condition-dependent parameters
    % duration of stimulus sequence in seconds
    h.Settings.totdur = 300; 
    % duration of trial in seconds
    h.Settings.trialdur = 0; % if 0, consecutive stimuli will occur with no gap
    % duration of stimulus in seconds
    h.Settings.stimdur = 0.8; % modified by oddball settings
    % Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.patternmethod = '';
    h.Settings.patternvalue = []; % one per stimdur
    % 'rand' or 'reg' spacing?
    h.Settings.stimdurtype = 'reg'; % not needed unless 'rand'
    % sampling rate
    h.Settings.fs = 96000; % don't change this
    % Binarual beats frequency: creates right ear frequency of f0+df
    h.Settings.df = 0; % 10Hz = alpha. Other options: 1Hz, 25Hz, 40Hz.
    % Monaural beats instead? 
    h.Settings.monaural = 0; 
    % attenuation level in decibels
    h.Settings.atten = -30; 
    % pitch
    h.Settings.f0 = 200; % Left ear carrier frequency (pitch)
    %intensity
    h.Settings.inten = 1; % value between 0 and 1

    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, pitch, channel
    h.Settings.oddballmethod = 'pitch'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.oddballvalue = [200 400];
    h.Settings.oddballtype = 'roving'; % options: 'roving', 'classical'
  
    
    %% SEQUENCE
    % Change probablity (CP): each condition is in rows
    h.Settings.oddprob = [
        % standard vs oddball
        0.8 0.2
        ];
    
    % index of row of oddprob that are standards and oddballs
    h.Settings.standardind = 1; % does not apply to 'roving oddball' design
    h.Settings.oddind = 2; % does not apply to 'roving oddball' design
    % keep oddball trials apart by at least sep_odd standards
    h.Settings.sep_odd = 2;
    % for sep_odd, which indices of h.Settings.oddballvalue to consider
    % each time? (each list will be separated separately)
    h.Settings.sep_odd_ind = {[1 2]};
    % for each set, ensure a number of leading standards 
    h.Settings.std_lead = 0;
    % number of minimum sets to randomised together
    h.Settings.n_set = 1; % 1 = use min set size
    % min number of oddballs within each CP condition
    h.Settings.n_odd = [1]; % overrides h.Settings.totdur
    % min number of oddballs per randomised set, per CP
    h.Settings.n_odd_set = [1]; % overrides h.Settings.totdur
    
    %% RESPONSE PARAMETERS
    % record responses during experiment? 0 or 1
    h.Settings.record_response = 1;
    % buttonpress options: key: keyboard inputs. Blank for no button press
    h.Settings.buttontype='key';
    % range of keyboard presses indicating a recordable response
    h.Settings.buttonopt = {'LeftArrow','RightArrow'}; 
    
    %% ADAPTIVE
    % adaptive staircase: meanings of the buttonopt
    h.Settings.adaptive.buttonfun = {'LeftArrow','RightArrow'}; 
    % adaptive staircase: corresponding signal values that would signify a
    % correct answer
    h.Settings.adaptive.signalval = [1 2];
    % starting level of adaptive staircase
    h.Settings.adaptive.startinglevel = max(h.Settings.oddballvalue)-min(h.Settings.oddballvalue); % for intensity, in dB (e.g. 10); for pitch, in Hz (e.g. 100)
    % adapt to omissions of response (not suitable for 2AFC tasks!!)
    h.Settings.adaptive.omit = 1; % 1 = omission is incorrect; 2 = omission is correct
    % which trials (or oddballs if oddonly selected) to start adaptive procedure if there is an omission?
    h.Settings.adaptive.startomit = 2;
    % adapt on every trial or only just before an oddball?
    h.Settings.adaptive.oddonly = 1;
    % max number of trials after oddball that subject must respond (otherwise counts as omitted response)
    h.Settings.adaptive.resptrials = 4;
    
    % number of trials of plot
    h.Settings.plottrials=0;
    

end

