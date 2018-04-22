function h = ALT(h,opt)

switch opt
    
    case 'setoptions'
        
    % settings options
    h.SettingsOptions = {'Ascend','Adaptive','ALT_classical','ALT_assoc'};
    
    case 'Ascend'

    % set general options
    h = setgeneral(h);
    
    % FILENAME OF SEQUENCE CREATION FUNCTION (without .m)
    h.SeqFun = 'SimpleSequence';
    
    %% TRIALS or CONTINUOUS?
    h.Settings.design = 'trials';
    % if continuous, how many trials ahead should be in the player schedule?
    % (applied to stimulation via soundcard only)
    h.Settings.ntrialsahead = 0;  %0 = all trials
    
    %% EXPERIMENTAL CONDIITIONS
    % name the settings that define orthogonal condtions at a different row
    h.Settings.conds = {};
    
    %% Output options
    % save sinwave from all trials as part of stim sequence file
    %h.Settings.savesinwave = 0;
    
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
    h.Settings.pauseeachblock = 0; % pause after every block?
    % names of any audiofiles
    h.Settings.audiofile = {};
    
    %% Condition-independent stimulus parameters - can be superceded by condition-dependent parameters
    % duration of stimulus sequence in seconds
    h.Settings.totdur = 0; 
    % duration of trial in seconds
    h.Settings.trialdur = 2; % if 0, consecutive stimuli will occur with no gap
    % Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.patternmethod = '';
    h.Settings.patternvalue = []; % one per stimdur
    % 'rand' or 'reg' spacing?
    h.Settings.stimdurtype = 'reg'; % not needed unless 'rand'
    % Binarual beats frequency: creates right ear frequency of f0+df
    %h.Settings.df = 0; % 10Hz = alpha. Other options: 1Hz, 25Hz, 40Hz.
    % Monaural beats instead? 
    %h.Settings.monaural = 0; 
    % attenuation level in decibels
    %h.Settings.atten = -30; 
    % pitch
    %h.Settings.f0 = 200; % Left ear carrier frequency (pitch)
    %intensity
    h.Settings.inten = 0; % value in dB
    h.Settings.maxinten = 0; % value in dB
    %intensity
    h.Settings.inten_diff = []; % value between 2 and 1000mA for Digitimer DS8R
    % Tactile: number of pulses per trial
    %h.Settings.nstim_trial = 1; % set to zero to be determined by stimdur
    % Tactile: within-trial frequency (Hz) 
    %h.Settings.t_freq=[]; 
    
   
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, pitch, channel
    h.Settings.oddballmethod = ''; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.oddballvalue = [];
    h.Settings.oddballtype = ''; % options: 'roving', 'classical'

    
    %% SEQUENCE
    h.Settings.ntrials = 500;
    % Change probablity (CP): each condition is in rows
    %h.Settings.oddprob = [
        % standard (left) vs oddball (right)=
    %    ];
    % index of row of oddprob that are standards and oddballs
    %h.Settings.standardind = 1; % does not apply to 'roving oddball' design
    % keep oddball trials apart by at least sep_odd standards
    %h.Settings.sep_odd = [2 1]; % for each CP condition
    % for each set, ensure a number of leading standards 
    %h.Settings.std_lead = [0 0]; % for each CP condition
    % number of minimum sets to randomise together
    %h.Settings.n_set = [1 1]; % 1 = use min set size, for each CP condition
    % min number of oddballs within each CP condition
    %h.Settings.n_odd = [20, 20]; % overrides h.Settings.totdur
    % min number of oddballs per randomised set, per CP
    %h.Settings.n_odd_set = [5, 5]; % overrides h.Settings.totdur
    
    %% RESPONSE PARAMETERS
    % record responses during experiment? 0 or 1
    h.Settings.record_response = 1;
    % buttonpress options: key: keyboard inputs. Blank for no button press
    h.Settings.buttontype='key';
    % range of keyboard presses indicating a recordable response
    h.Settings.buttonopt = {'UpArrow','DownArrow'}; 
    
    %% THRESHOLDING
    % starting level and step size
    h.Settings.threshold.type = 'intensity'; % for intensity
    h.Settings.threshold.startinglevel = h.Settings.atten-100; % for intensity, in dB (e.g. 10); for pitch, in Hz (e.g. 100) 
    h.Settings.threshold.step = 5;
    h.Settings.threshold.signalval = [1 2]; % 1 = carrying on increasing; 2 = decrease
    h.Settings.threshold.maxinten = 0; % 0 is the max
    
case 'Adaptive'

    % set general options
    h = setgeneral(h);
    
    % FILENAME OF SEQUENCE CREATION FUNCTION (without .m)
    h.SeqFun = 'CreateSequence';
    
    %% TRIALS or CONTINUOUS?
    h.Settings.design = 'trials';
    % if continuous, how many trials ahead should be in the player schedule?
    % (applied to stimulation via soundcard only)
    h.Settings.ntrialsahead = 0;  %0 = all trials
    
    %% EXPERIMENTAL CONDIITIONS
    % name the settings that define orthogonal condtions at a different row
    h.Settings.conds = {};
    
    %% Output options
    % save sinwave from all trials as part of stim sequence file
    %h.Settings.savesinwave = 0;
    
    %% BLOCKING/RUN OPTIONS
    % 'divide' = equally divide trials by nblocks; 
    % 'cond' = separate block for each condition
    h.Settings.blockopt = 'divide';
    % further options for 'divide':
        % number of blocks (containing multiple conditions)
        %h.Settings.nblocks = 1; % must integer-divide each value in h.Settings.cond_rep_init
        %distribute conditions equally among blocks
    %    h.Settings.distblocks = 1;
    % options to start sequence at beginning of every run
    % 'msgbox', 'labjack', 'buttonpress', 'audio' - can have more than one in
    % cell array
    h.Settings.blockstart = {'buttonpress'}; % audio,labjack,audio
    h.Settings.pauseeachblock = 0; % pause after every block?
    % names of any audiofiles
    h.Settings.audiofile = {};
    
    %% Condition-independent stimulus parameters - can be superceded by condition-dependent parameters
    % duration of stimulus sequence in seconds
    h.Settings.totdur = 0; 
    % duration of trial in seconds
    h.Settings.trialdur = 3; % if 0, consecutive stimuli will occur with no gap
    % Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.nstim_trial = 2; % set to zero to be determined by stimdur
    % Tactile: within-trial frequency (Hz) 
    h.Settings.wait=[1 0]; % one value per nstim 
    
    %% first stimulus: audio
    % duration of stimulus in seconds
    h.Settings.stim(1).dur = 0.3; % modified by oddball settings
    h.Settings.stim(1).patternmethod = '';% Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.stim(1).patternvalue = {}; % one per stimdur in each cell; one cell per oddball value
    h.Settings.stim(1).durtype = 'reg'; % not needed unless 'rand'
    h.Settings.stim(1).inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).f0 = 500; % pitch
    h.Settings.stim(1).inten_type = 'dB'; % either 'dB' or 'abs'
    h.Settings.stim(1).df = 0;
    h.Settings.stim(1).atten = {'inten_diff',1}; % attenuation level in decibels OR cell with variable and multiplier
    h.Settings.stim(1).attenchan = [1 2]; % apply attenuation (e.g. during thresholding) to these chans
    h.Settings.stim(1).control='PsychPortAudio'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    h.Settings.stim(1).chan = [1 2]; 
    h.Settings.stim(1).nrchannels = 2; % total number of channels, e.g. on sound card
    h.Settings.stim(1).Tukey = 0.25; % Apply Tukey window?
    h.Settings.stim(1).Tukeytype = 2; % 1 = apply to each tone within pattern; 2 = apply to whole pattern

    %% second stimulus: audio
    % duration of stimulus in seconds
    h.Settings.stim(2).dur = 0.3; % modified by oddball settings
    h.Settings.stim(2).patternmethod = '';% Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.stim(2).patternvalue = []; % one per stimdur. Can leave empty if using "index" as oddball method with GUI to set inten_diff
    h.Settings.stim(2).durtype = 'reg'; % not needed unless 'rand'
    h.Settings.stim(2).inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(2).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(2).maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(2).f0 = 500; % pitch
    h.Settings.stim(2).inten_type = 'dB'; % either 'dB' or 'abs'
    h.Settings.stim(2).df = 0;
    h.Settings.stim(2).atten = {'inten_diff',1}; % attenuation level in decibels
    h.Settings.stim(2).attenchan = [1 2]; % apply attenuation (e.g. during thresholding) to these chans
    h.Settings.stim(2).control='PsychPortAudio'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    h.Settings.stim(2).chan = [1 2]; 
    h.Settings.stim(2).nrchannels = 2; % total number of channels, e.g. on sound card
    h.Settings.stim(2).Tukey = 0.25; % Apply Tukey window?
    h.Settings.stim(2).Tukeytype = 2; % 1 = apply to each tone within pattern; 2 = apply to whole pattern
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, index, pitch, channel
    h.Settings.oddballmethod = 'intensity'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.oddballvalue = {[1 2]}; % values to go into h.Seq.signal. One per oddprob row, or leave blank if determined from GUI
    h.Settings.oddballtype = 'classical'; % options: 'roving', 'classical'

    %% SEQUENCE
    % Change probablity (CP): each condition is in rows
    h.Settings.oddprob = [
        % standard (left) vs oddball (right)
        0.5 0.5
        ];
    % index of row of oddprob that are standards and oddballs. Can be
    % overridden by h.Settings.oddballvalue if using index
    h.Settings.standardind = 1; 
    h.Settings.oddind = 2; 
    % keep oddball trials apart by at least sep_odd standards
    h.Settings.sep_odd = [0]; % for each CP condition
    % for sep_odd, which indices of h.Settings.oddballvalue to consider
    % each time? (each list will be considered separately)
    h.Settings.sep_odd_ind = {[1 2]};
    % for each set, ensure a number of leading standards 
    h.Settings.std_lead = [0]; % for each CP condition
    % number of sets to randomise together
    h.Settings.n_set = []; % Leave blank to calculate automatically; or one nunmber per CP condition
    % min number of oddballs within each CP condition
    h.Settings.n_odd = [100*12]; % overrides h.Settings.totdur
    % min number of oddballs per randomised set, per CP
    h.Settings.n_odd_set = [100*12]; % overrides h.Settings.totdur
    % randomise sets?
    h.Settings.rand_set = [0]; 
    
    %% RESPONSE PARAMETERS
    % record responses during experiment? 0 or 1
    h.Settings.record_response = 1;
    % how to record responses?
    h.Settings.record_response_type = {'all'}; %options: 'all','thistrial','previoustrial'
    %intensity difference
    % buttonpress options: key: keyboard inputs. Blank for no button press
    h.Settings.buttontype='key';
    % range of keyboard presses indicating a recordable response
    h.Settings.buttonopt = {'DownArrow','UpArrow'}; 
    
    
    %% ADAPTIVE: General
    % which ones to run? (i.e. indices of h.Settings.adaptive)
    h.Settings.adaptive_general.adapttypes = [1];
    % alternate or randomise runs over types? Alt must have equal number of
    % runs for each adapttype
    h.Settings.adaptive_general.seqtype = 'rand'; % 'alt' or 'rand'
    h.Settings.adaptive_general.seqrandblocksize = 12; % should divide the number of trials in a set
    h.Settings.adaptive_general.selectcond.cp = [1]; % which CP condition to run adaptive on?
    h.Settings.adaptive_general.stim = 2; % which stim to run adaptive on?
    
    %% ADAPTIVE 1
    h.Settings.adaptive(1).type = 'discrim';
    h.Settings.adaptive(1).updown = [1 2];
    % how many of each to run?
    h.Settings.adaptive(1).nRuns = 12*100;
    % max number of thresh estimates to average over to get overall estimate
    h.Settings.adaptive(1).av_thresh = [50,75,100];
    % number of trials each run
    h.Settings.adaptive(1).trialsperrun = 1;
    % adaptive staircase: meanings of the buttonopt
    h.Settings.adaptive(1).buttonfun = {'DownArrow','UpArrow'}; 
    % adaptive staircase: corresponding signal values that would signify a
    % correct answer
    h.Settings.adaptive(1).signalval = [2 1];
    % reversals
    h.Settings.adaptive(1).reversals = [4;8;12];
    % stepsize
    h.Settings.adaptive(1).stepsize = [2;sqrt(2);sqrt(sqrt(2))];
    % steptype 0 = multiple/divide by stepsize; steptype 1 = add/subtract
    h.Settings.adaptive(1).steptype = 0;
    % stepdir -1 = level decreases intensity; stepdir 1 = level increases intensity
    h.Settings.adaptive(1).stepdir = -1;
    % starting level of adaptive staircase
    h.Settings.adaptive(1).startinglevel = -5; % should be a DIFFERENCE value in dB. Keep small as it will increase naturally over time.
    % adapt to omissions of response (not suitable for 2AFC tasks, so set to 0)
    h.Settings.adaptive(1).omit = 0; % 1 = omission is incorrect; 2 = omission is correct
    % which trials (or oddballs if oddonly selected) to start adaptive procedure if there is an omission?
    h.Settings.adaptive(1).startomit = 0;
    % adapt on every trial or only just before an oddball?
    %h.Settings.adaptive.oddonly = 1;
    % max number of trials after oddball that subject must respond (otherwise counts as omitted response)
    %h.Settings.adaptive.resptrials = 4;
    % number of reversals to average over to calculate threshold.
    h.Settings.adaptive(1).reversalForthresh = 6;
    % use mean from the first X responses of each type (high and low)
    %h.Settings.adaptive(1).getmeanfromresponses = 6;
    % maximum amount to adjust the mean if their responses are very
    % incorrect (should be a small fraction, e.g. 1/5th, of the stimulus intensity)
    %h.Settings.adaptive(1).meanadjustmax = 10;
    % maximum amount of the difference value (should be a small fraction, e.g. 1/5th, of the stimulus intensity)
    h.Settings.adaptive(1).levelmax = -50; % should be a DIFFERENCE value in dB. 
    
    
    case 'ALT_classical'

    % set general options
    h = setgeneral(h);
    
    % FILENAME OF SEQUENCE CREATION FUNCTION (without .m)
    h.SeqFun = 'CreateSequence';
    
    %% TRIALS or CONTINUOUS?
    h.Settings.design = 'trials';
    % if continuous, how many trials ahead should be in the player schedule?
    % (applied to stimulation via soundcard only)
    h.Settings.ntrialsahead = 0;  %0 = all trials
    
    %% EXPERIMENTAL CONDIITIONS
    % name the settings that define orthogonal condtions at a different row
    h.Settings.conds = {'oddprob'};
    
    %% Output options
    % save sinwave from all trials as part of stim sequence file
    %h.Settings.savesinwave = 0;
    
    %% BLOCKING/RUN OPTIONS
    % 'divide' = equally divide trials by nblocks; 
    % 'cond' = separate block for each condition
    h.Settings.blockopt = 'divide';
    % further options for 'divide':
        % number of blocks (containing multiple conditions)
    %    h.Settings.nblocks = 2; % must integer-divide each value in h.Settings.cond_rep_init
        %distribute conditions equally among blocks
    %    h.Settings.distblocks = 1;
    % options to start sequence at beginning of every run
    % 'msgbox', 'labjack', 'buttonpress', 'audio' - can have more than one in
    % cell array
    h.Settings.blockstart = {'buttonpress'}; % audio,labjack,audio
    h.Settings.pauseeachblock = 0; % pause after every block?
    % names of any audiofiles
    h.Settings.audiofile = {};
    
    %% Condition-independent stimulus parameters - can be superceded by condition-dependent parameters
    % duration of stimulus sequence in seconds
    h.Settings.totdur = 0; 
    % duration of trial in seconds
    h.Settings.trialdur = 2; % if 0, consecutive stimuli will occur with no gap
    % Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.patternmethod = '';
    h.Settings.patternvalue = []; % one per stimdur
    % 'rand' or 'reg' spacing?
    h.Settings.stimdurtype = 'reg'; % not needed unless 'rand'
    % Binarual beats frequency: creates right ear frequency of f0+df
    %h.Settings.df = 0; % 10Hz = alpha. Other options: 1Hz, 25Hz, 40Hz.
    % Monaural beats instead? 
    %h.Settings.monaural = 0; 
    % attenuation level in decibels
    %h.Settings.atten = -30; 
    % pitch
    %h.Settings.f0 = 200; % Left ear carrier frequency (pitch)
    %intensity
    h.Settings.inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    % Tactile: number of pulses per trial
    %h.Settings.nstim_trial = 1; % set to zero to be determined by stimdur
    % Tactile: within-trial frequency (Hz) 
    %h.Settings.t_freq=[]; 
    
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, index, pitch, channel
    h.Settings.oddballmethod = 'index'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.oddballvalue = {[1 2], [1 2], [2 1]}; % values to go into h.Seq.signal. One per oddprob row, or leave blank if determined from GUI
    h.Settings.oddballtype = 'classical'; % options: 'roving', 'classical'

    %% SEQUENCE
    % Change probablity (CP): each condition is in rows
    h.Settings.oddprob = [
        % standard (left) vs oddball (right)
        0.8 0.2
        0.5 0.5
        0.8 0.2
        ];
    % index of row of oddprob that are standards and oddballs. Can be
    % overridden by h.Settings.oddballvalue if using index
    h.Settings.standardind = 1; 
    h.Settings.oddind = 2; 
    % keep oddball trials apart by at least sep_odd standards
    h.Settings.sep_odd = [2 0 2]; % for each CP condition
    % for sep_odd, which indices of h.Settings.oddballvalue to consider
    % each time? (each list will be considered separately)
    h.Settings.sep_odd_ind = {[1 2],[1 2],[1 2]};
    % for each set, ensure a number of leading standards 
    h.Settings.std_lead = [0 0 0]; % for each CP condition
    % number of sets to randomise together
    h.Settings.n_set = []; % Leave blank to calculate automatically; or one nunmber per CP condition
    % min number of oddballs within each CP condition
    h.Settings.n_odd = [12, 30, 12]; % overrides h.Settings.totdur
    % min number of oddballs per randomised set, per CP
    h.Settings.n_odd_set = [4, 10, 4]; % overrides h.Settings.totdur
    % randomise sets?
    h.Settings.rand_set = [1 1 1]; 
    
    %% RESPONSE PARAMETERS
    % record responses during experiment? 0 or 1
    h.Settings.record_response = 1;
    % how to record responses?
    h.Settings.record_response_type = {'all'}; %options: 'all','thistrial','previoustrial'
    %intensity difference
    % buttonpress options: key: keyboard inputs. Blank for no button press
    h.Settings.buttontype='key';
    % range of keyboard presses indicating a recordable response
    h.Settings.buttonopt = {'DownArrow','UpArrow'}; 
    
    %% THRESHOLDING
    % starting level and step size
    %h.Settings.threshold.startinglevel = 2; % for intensity)
    %h.Settings.threshold.step = 2;
    
    
    case 'ALT_assoc'

    % set general options
    h = setgeneral(h);
    
    % FILENAME OF SEQUENCE CREATION FUNCTION (without .m)
    h.SeqFun = 'CreateSequence';
    
    %% TRIALS or CONTINUOUS?
    h.Settings.design = 'trials';
    % if continuous, how many trials ahead should be in the player schedule?
    % (applied to stimulation via soundcard only)
    h.Settings.ntrialsahead = 0;  %0 = all trials
    
    %% EXPERIMENTAL CONDIITIONS
    % name the settings that define orthogonal conditions at a different row
    h.Settings.conds = {'oddprob'};
    
    %% Output options
    % save sinwave from all trials as part of stim sequence file
    %h.Settings.savesinwave = 0;
    
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
    h.Settings.pauseeachblock = 0; % pause after every block?
    % names of any audiofiles
    h.Settings.audiofile = {};
    
    %% Condition-independent stimulus parameters - can be superceded by condition-dependent parameters
    % duration of stimulus sequence in seconds
    h.Settings.totdur = 0; 
    % duration of trial in seconds
    h.Settings.trialdur = 3; % if 0, consecutive stimuli will occur with no gap
    % Tactile: number of pulses per trial
    h.Settings.nstim_trial = 2; % set to zero to be determined by stimdur
    % Tactile: within-trial frequency (Hz) 
    h.Settings.wait=[1 0]; % one value per nstim 
    
    %% first stimulus: audio
    % duration of stimulus in seconds
    h.Settings.stim(1).dur = [0.15 0.15]; % modified by oddball settings
    h.Settings.stim(1).patternmethod = 'pitch';% Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.stim(1).patternvalue = {[540 500],[460 500]}; % one per stimdur in each cell; one cell per oddball value
    h.Settings.stim(1).durtype = 'reg'; % not needed unless 'rand'
    h.Settings.stim(1).inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).f0 = 500; % pitch
    h.Settings.stim(1).inten_type = 'dB'; % either 'dB' or 'abs'
    h.Settings.stim(1).df = 0;
    h.Settings.stim(1).atten = {'inten_diff',1}; % attenuation level in decibels
    h.Settings.stim(1).attenchan = [1 2]; % apply attenuation (e.g. during thresholding) to these chans
    h.Settings.stim(1).control='PsychPortAudio'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    h.Settings.stim(1).chan = [1 2]; 
    h.Settings.stim(1).nrchannels = 2; % total number of channels, e.g. on sound card
    h.Settings.stim(1).Tukey = 0.25; % Apply Tukey window?
    h.Settings.stim(1).Tukeytype = 2; % 1 = apply to each tone within pattern; 2 = apply to whole pattern

    %% second stimulus: audio
    % duration of stimulus in seconds
    h.Settings.stim(2).dur = 0.3; % modified by oddball settings
    h.Settings.stim(2).patternmethod = 'intensity';% Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.stim(2).patternvalue = []; % one per stimdur. Can leave empty if using "index" as oddball method with GUI to set inten_diff
    h.Settings.stim(2).durtype = 'reg'; % not needed unless 'rand'
    h.Settings.stim(2).inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(2).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(2).maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(2).f0 = 500; % pitch
    h.Settings.stim(2).inten_type = 'dB'; % either 'dB' or 'abs'
    h.Settings.stim(2).df = 0;
    h.Settings.stim(2).atten = {'inten_diff',1}; % attenuation level in decibels
    h.Settings.stim(2).attenchan = [1 2]; % apply attenuation (e.g. during thresholding) to these chans
    h.Settings.stim(2).control='PsychPortAudio'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    h.Settings.stim(2).chan = [1 2]; 
    h.Settings.stim(2).nrchannels = 2; % total number of channels, e.g. on sound card
    h.Settings.stim(2).Tukey = 0.25; % Apply Tukey window?
    h.Settings.stim(2).Tukeytype = 2; % 1 = apply to each tone within pattern; 2 = apply to whole pattern
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, index, pitch, channel
    h.Settings.oddballmethod = 'index'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.oddballvalue = {[1 2], [1 2],  [1 2], [2 1], [1 2],  [1 2], [1 2], [2 1], [1 2],  [1 2], [1 2], [2 1]}; % values to go into h.Seq.signal. One per oddprob row, or leave blank if determined from GUI
    h.Settings.oddballtype = 'classical'; % options: 'roving', 'classical'

    %% SEQUENCE
    % Change probablity (CP): each condition is in rows
    h.Settings.oddprob = [
        % standard (left) vs oddball (right)
        0.8 0.2
        0.5 0.5
        0.8 0.2
        0.8 0.2
        0.5 0.5
        0.8 0.2
        0.8 0.2
        0.5 0.5
        0.8 0.2
        ];
    % index of row of oddprob that are standards and oddballs. Can be
    % overridden by h.Settings.oddballvalue if using index
    h.Settings.standardind = 1; 
    h.Settings.oddind = 2; 
    % keep oddball trials apart by at least sep_odd standards
    h.Settings.sep_odd = [2 0 2 2 0 2 2 0 2]; % for each CP condition
    % for sep_odd, which indices of h.Settings.oddballvalue to consider
    % each time? (each list will be considered separately)
    h.Settings.sep_odd_ind = {[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2]};
    % for each set, ensure a number of leading standards 
    h.Settings.std_lead = [0 0 0 0 0 0 0 0 0]; % for each CP condition
    % number of sets to randomise together
    h.Settings.n_set = []; % Leave blank to calculate automatically; or one nunmber per CP condition
    % min number of oddballs within each CP condition
    h.Settings.n_odd = [4, 4, 4,  8, 8, 8,  12, 12, 12]; % overrides h.Settings.totdur
    % min number of oddballs per randomised set, per CP
    h.Settings.n_odd_set = h.Settings.n_odd; % overrides h.Settings.totdur
    % randomise sets?
    h.Settings.rand_set = [1 1 1 1 1 1 1 1 1]; 
    % condition numbers
    h.Settings.condnum = [
        1 2
        3 4
        5 6
        1 2
        3 4
        5 6
        1 2
        3 4
        5 6
        ]; 
    
    %% associative learning experiments
    h.Settings.assoc.pairing = {
      % cue 1  cue 2
        [1 3], [2 4]; %pairing option 1 (e.g. most porbable in some blocks)
        [2 4], [1 3]; %pairing option 2 (e.g. least probable)
        }; % more than one number per pairing/cue: will be randomly assigned to each cue
    % if there is more than on number per pairing/cue, what is their
    % percent representation?
    h.Settings.assoc.probstim = [0.5 0.5];
    % for each condnum, which pair to use?
    h.Settings.assoc.pair = [1 2 1 2 2 1];
    % for each stimtype (unique value) within h.Settings.assoc.pairing, 
    % what inten_diff multiplier to use?
    h.Settings.assoc.intenstim = [1 -1 2 -2];
    
    %% RESPONSE PARAMETERS
    % record responses during experiment? 0 or 1
    h.Settings.record_response = 1;
    % how to record responses?
    h.Settings.record_response_type = {'all'}; %options: 'all','thistrial','previoustrial'
    %intensity difference
    % buttonpress options: key: keyboard inputs. Blank for no button press
    h.Settings.buttontype='key';
    % range of keyboard presses indicating a recordable response
    h.Settings.buttonopt = {'DownArrow','UpArrow'}; 
    
    %% THRESHOLDING
    % starting level and step size
    %h.Settings.threshold.startinglevel = 2; % for intensity)
    %h.Settings.threshold.step = 2;
    
    
end

function h = setgeneral(h)

%% EQUIPMENT CONTROL
% record EEG, NS: netstation, BV: brainvision, 'serial': serial port
% serial port
%h.Settings.serial = 'COM1';
%h.Settings.record_EEG='serial';
%h.Settings.EEGport = 'COM3'; % only needed for 'serial' EEG triggers
%h.Settings.EEGMarkPattern = 1; % mark EEG for every change in stimulus pattern (0 = start of trial only)
h.Settings.labjack=0; % Use labjack for controlling any equipment?
%h.Settings.stimcontrol='PsychPortAudio'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
%h.Settings.labjack_DACport = 0;
%h.Settings.DAC_multiply = 0; % multiply DAC output by this (e.g. to get mA on DS8R)
%h.Settings.stimchan = [1 2]; 
%h.Settings.stimchanforLJ = 0;
%h.Settings.nrchannels = 2; % total number of channels, e.g. on sound card
% channels on stimulator to use; use differenr rows for different pairs
% (e.g. for different conditions). If labjack, this can refer to output
% port(s).

% Tactile: number of pulses in a train
%h.Settings.npulses_train = 1; % set to zero to be determined by stimdur
% Tactile: within-train frequency (Hz)
%h.Settings.p_freq=[];
% Use timer to control frequency of labjack outputs? Otherwise uses software timing.
%h.Settings.labjack_timer=1; % must use if there are a large number of pulses (e.g. >4)

%% AUDIO

% duration of stimulus in seconds
%h.Settings.stimdur = 0.3; % modified by oddball settings
% pitch
%h.Settings.f0 = 500; 
%intensity
%h.Settings.inten = 0; % value between 0 and 1, or if decibels value of <=0
%h.Settings.inten_type = 'dB'; % either 'dB' or 'abs'
%h.Settings.df = 0;
% Audio sampling rate 
h.Settings.fs = 96000; % don't change this
% attenuation level in decibels
%h.Settings.atten = 0; 
%h.Settings.attenchan = [1 2]; % apply attenuation (e.g. during thresholding) to these chans
% stimulate responses?
h.Settings.simulate_response = 0;