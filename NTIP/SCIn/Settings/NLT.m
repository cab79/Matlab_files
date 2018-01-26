function h = NLT(h,opt)

switch opt
    
    case 'setoptions'
        
    % settings options
    h.SettingsOptions = {'Ascend','Adapt_eachtrial','Adapt_overtrials','NLT_roving','NLT_classical','Complex','NLT_complex','NLT_complex_entrain'};
    
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
    % names of any audiofiles
    h.Settings.audiofile = {};
    
    %% Condition-independent stimulus parameters - can be superceded by condition-dependent parameters
    % duration of stimulus sequence in seconds
    h.Settings.totdur = 0; 
    % duration of trial in seconds
    h.Settings.trialdur = 2; % if 0, consecutive stimuli will occur with no gap
    % duration of stimulus in seconds
    h.Settings.stimdur = 0; % modified by oddball settings
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
    h.Settings.inten = []; % value between 2 and 1000mA for Digitimer DS8R
    % Tactile: number of pulses per trial
    h.Settings.nstim_trial = 1; % set to zero to be determined by stimdur
    % Tactile: within-trial frequency (Hz) 
    h.Settings.t_freq=[]; 
    
   
    
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
    h.Settings.buttonopt = {'DownArrow'}; 
    
    %% THRESHOLDING
    % starting level and step size
    h.Settings.threshold.startinglevel = 0; % for intensity)
    h.Settings.threshold.step = 2;
    
case 'Adapt_eachtrial'

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
    h.Settings.audiofile = {};
    
    %% Condition-independent stimulus parameters - can be superceded by condition-dependent parameters
    % duration of stimulus sequence in seconds
    h.Settings.totdur = 0; 
    % duration of trial in seconds
    h.Settings.trialdur = 5; % if 0, consecutive stimuli will occur with no gap
    % duration of stimulus in seconds
    h.Settings.stimdur = 0; % modified by oddball settings
    % Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.patternmethod = '';
    h.Settings.patternvalue = []; % leave empty to product values from baseintensity and adaptive
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
    h.Settings.inten = []; % value between 2 and 1000mA for Digitimer DS8R
    % Tactile: number of pulses per trial
    h.Settings.nstim_trial = 2; % set to zero to be determined by stimdur
    % Tactile: within-trial frequency (Hz) 
    h.Settings.t_freq=1; 
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, pitch, channel
    h.Settings.oddballmethod = 'intensity'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.oddballvalue = []; % leave blank if determined from GUI
    h.Settings.oddballtype = 'classical'; % options: 'roving', 'classical'

    %% SEQUENCE
    % Change probablity (CP): each condition is in rows
    h.Settings.oddprob = [
        % standard (left) vs oddball (right)
        0.5 0.5
        ];
    % index of row of oddprob that are standards and oddballs
    h.Settings.standardind = 1; 
    h.Settings.oddind = 2; 
    % keep oddball trials apart by at least sep_odd standards
    h.Settings.sep_odd = [0]; % for each CP condition
    % for sep_odd, which indices of h.Settings.oddballvalue to consider
    % each time? (each list will be considered separately)
    h.Settings.sep_odd_ind = {};
    % for each set, ensure a number of leading standards 
    h.Settings.std_lead = [0]; % for each CP condition
    % min number of oddballs within each CP condition
    h.Settings.n_odd = [30]; % overrides h.Settings.totdur
    % number of minimum sets to randomise together
    h.Settings.n_set = [1]; % 1 = use min set size, for each CP condition
    % min number of oddballs per randomised set, per CP and per oddball
    % type
    h.Settings.n_odd_set = [30]; % overrides h.Settings.totdur
    % randomise sets?
    h.Settings.rand_set = 0; 
    
    %% RESPONSE PARAMETERS
    % record responses during experiment? 0 or 1
    h.Settings.record_response = 1;
    % how to record responses?
    h.Settings.record_response_type = {'thistrial','previoustrial'}; %options: 'all','thistrial','previoustrial'
    % buttonpress options: key: keyboard inputs. Blank for no button press
    h.Settings.buttontype='key';
    % range of keyboard presses indicating a recordable response
    h.Settings.buttonopt = {'LeftArrow','RightArrow'}; 
    
    %% ADAPTIVE
    % adaptive staircase: meanings of the buttonopt
    h.Settings.adaptive.buttonfun = {'LeftArrow','RightArrow'}; 
    % adaptive staircase: corresponding signal values that would signify a
    % correct answer
    h.Settings.adaptive.signalval = [2 1];
    % steptype 0 = multiple/divide by stepsize; steptype 1 = add/subtract
    h.Settings.steptype = 0;
    % stepdir -1 = level decreases intensity; stepdir 1 = level increases intensity
    h.Settings.stepdir = -1;
    % starting level of adaptive staircase
    h.Settings.adaptive.startinglevel = 10; % should be a DIFFERENCE value; for intensity, in dB (e.g. 10); for pitch, in Hz (e.g. 100); for tactile in mA
    % adapt to omissions of response (not suitable for 2AFC tasks, so set to 0)
    h.Settings.adaptive.omit = 0; % 1 = omission is incorrect; 2 = omission is correct
    % which trials (or oddballs if oddonly selected) to start adaptive procedure if there is an omission?
    h.Settings.adaptive.startomit = 0;
    % adapt on every trial or only just before an oddball?
    %h.Settings.adaptive.oddonly = 1;
    % max number of trials after oddball that subject must respond (otherwise counts as omitted response)
    %h.Settings.adaptive.resptrials = 4;
    % number of reversals to average over to calculate threshold.
    h.Settings.reversalForthresh = 6;
    
case 'Adapt_overtrials'

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
    h.Settings.audiofile = {};
    
    %% Condition-independent stimulus parameters - can be superceded by condition-dependent parameters
    % duration of stimulus sequence in seconds
    h.Settings.totdur = 0; 
    % duration of trial in seconds
    h.Settings.trialdur = 3; % if 0, consecutive stimuli will occur with no gap
    % duration of stimulus in seconds
    h.Settings.stimdur = 0; % modified by oddball settings
    % Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.patternmethod = '';
    h.Settings.patternvalue = []; % leave empty to product values from baseintensity and adaptive
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
    h.Settings.inten = []; % value between 2 and 1000mA for Digitimer DS8R
    % Tactile: number of pulses per trial
    h.Settings.nstim_trial = 1; % set to zero to be determined by stimdur
    % Tactile: within-trial frequency (Hz) 
    h.Settings.t_freq=[]; 
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, pitch, channel
    h.Settings.oddballmethod = 'intensity'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.oddballvalue = []; % leave blank if determined from GUI
    h.Settings.oddballtype = 'classical'; % options: 'roving', 'classical'

    %% SEQUENCE
    % Change probablity (CP): each condition is in rows
    h.Settings.oddprob = [
        % standard (left) vs oddball (right)
        0.5 0.5
        ];
    % index of row of oddprob that are standards and oddballs
    h.Settings.standardind = 1; 
    h.Settings.oddind = 2; 
    % keep oddball trials apart by at least sep_odd standards
    h.Settings.sep_odd = [0]; % for each CP condition
    % for sep_odd, which indices of h.Settings.oddballvalue to consider
    % each time? (each list will be considered separately)
    h.Settings.sep_odd_ind = {};
    % for each set, ensure a number of leading standards 
    h.Settings.std_lead = [0]; % for each CP condition
    % min number of oddballs within each CP condition
    h.Settings.n_odd = [60]; % overrides h.Settings.totdur
    % number of minimum sets to randomise together
    h.Settings.n_set = [1]; % 1 = use min set size, for each CP condition
    % min number of oddballs per randomised set, per CP and per oddball
    % type
    h.Settings.n_odd_set = [60]; % overrides h.Settings.totdur
    % randomise sets?
    h.Settings.rand_set = 0; 
    
    %% RESPONSE PARAMETERS
    % record responses during experiment? 0 or 1
    h.Settings.record_response = 1;
    % how to record responses?
    h.Settings.record_response_type = {'thistrial','previoustrial'}; %options: 'all','thistrial','previoustrial'
    % buttonpress options: key: keyboard inputs. Blank for no button press
    h.Settings.buttontype='key';
    % range of keyboard presses indicating a recordable response
    h.Settings.buttonopt = {'LeftArrow','RightArrow'}; 
    
    %% ADAPTIVE: General
    % Order of running
    h.Settings.adaptive_general.order = [1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2];
    
    %% ADAPTIVE 1
    h.Settings.adaptive(1).type = 'detect';
    h.Settings.adaptive(1).updown = [1 1];
    % number of trials each run
    h.Settings.adaptive(1).trialsperrun = 12;
    % adaptive staircase: meanings of the buttonopt
    h.Settings.adaptive(1).buttonfun = {'LeftArrow','RightArrow'}; 
    % adaptive staircase: corresponding signal values that would signify a
    % correct answer
    h.Settings.adaptive(1).signalval = [1 2];
    % reversals
    h.Settings.adaptive(1).reversals = [4;8];
    % stepsize
    h.Settings.adaptive(1).stepsize = [2;sqrt(2)];
    % steptype 0 = multiple/divide by stepsize; steptype 1 = add/subtract
    h.Settings.adaptive(1).steptype = 0;
    % stepdir -1 = level decreases intensity; stepdir 1 = level increases intensity
    h.Settings.adaptive(1).stepdir = -1;
    % starting level of adaptive staircase
    h.Settings.adaptive(1).startinglevel = 2; % should be a DIFFERENCE value in mA. Keep small as it will increase naturally over time.
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
    h.Settings.adaptive(1).levelmax = 10; % should be a DIFFERENCE value in mA. Keep small as it will increase naturally over time.
    
    %% ADAPTIVE 2
    h.Settings.adaptive(2).type = 'discrim';
    h.Settings.adaptive(2).updown = [1 2];
    % number of trials each run
    h.Settings.adaptive(2).trialsperrun = 12;
    % adaptive staircase: meanings of the buttonopt
    h.Settings.adaptive(2).buttonfun = {'LeftArrow','RightArrow'}; 
    % adaptive staircase: corresponding signal values that would signify a
    % correct answer
    h.Settings.adaptive(2).signalval = [1 2];
    % reversals
    h.Settings.adaptive(2).reversals = [4;8];
    % stepsize
    h.Settings.adaptive(2).stepsize = [2;sqrt(2)];
    % steptype 0 = multiple/divide by stepsize; steptype 1 = add/subtract
    h.Settings.adaptive(2).steptype = 0;
    % stepdir -1 = level decreases intensity; stepdir 1 = level increases intensity
    h.Settings.adaptive(2).stepdir = -1;
    % starting level of adaptive staircase
    h.Settings.adaptive(2).startinglevel = 2; % should be a DIFFERENCE value in mA. Keep small as it will increase naturally over time.
    % adapt to omissions of response (not suitable for 2AFC tasks, so set to 0)
    h.Settings.adaptive(2).omit = 0; % 1 = omission is incorrect; 2 = omission is correct
    % which trials (or oddballs if oddonly selected) to start adaptive procedure if there is an omission?
    h.Settings.adaptive(2).startomit = 0;
    % adapt on every trial or only just before an oddball?
    %h.Settings.adaptive.oddonly = 1;
    % max number of trials after oddball that subject must respond (otherwise counts as omitted response)
    %h.Settings.adaptive.resptrials = 4;
    % number of reversals to average over to calculate threshold.
    h.Settings.adaptive(2).reversalForthresh = 6;
    % use mean from the first X responses of each type (high and low)
    %h.Settings.adaptive(2).getmeanfromresponses = 6;
    % maximum amount to adjust the mean if their responses are very
    % incorrect (should be a small fraction, e.g. 1/5th, of the stimulus intensity)
    %h.Settings.adaptive(2).meanadjustmax = 10;
    % maximum amount of the difference value (should be a small fraction, e.g. 1/5th, of the stimulus intensity)
    h.Settings.adaptive(2).levelmax = 10; % should be a DIFFERENCE value in mA. Keep small as it will increase naturally over time.
    
    case 'NLT_roving'

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
    % names of any audiofiles
    h.Settings.audiofile = {};
    
    %% Condition-independent stimulus parameters - can be superceded by condition-dependent parameters
    % duration of stimulus sequence in seconds
    h.Settings.totdur = 0; 
    % duration of trial in seconds
    h.Settings.trialdur = 0.4; % if 0, consecutive stimuli will occur with no gap
    % duration of stimulus in seconds
    h.Settings.stimdur = 0; % modified by oddball settings
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
    h.Settings.inten = []; % value between 2 and 1000mA for Digitimer DS8R. leave blank if determined from GUI, threshold or adaptive
    % Tactile: number of pulses per trial
    h.Settings.nstim_trial = 1; % set to zero to be determined by stimdur
    % Tactile: within-trial frequency (Hz) 
    h.Settings.t_freq=0; 
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, pitch, channel
    h.Settings.oddballmethod = 'intensity'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.oddballvalue = []; % leave blank if determined from GUI
    h.Settings.oddballtype = 'roving'; % options: 'roving', 'classical'

    %% SEQUENCE
    % Change probablity (CP): each condition is in rows
    h.Settings.oddprob = [
        % standard (left) vs oddball (right)
        0.9 0.1
        0.7 0.3
        0.5 0.5
        ];
    % index of row of oddprob that are standards and oddballs
    h.Settings.standardind = 1; 
    h.Settings.oddind = 2; 
    % keep oddball trials apart by at least sep_odd standards
    h.Settings.sep_odd = [2 2 1]; % for each CP condition
    % for sep_odd, which indices of h.Settings.oddballvalue to consider
    % each time? (each list will be considered separately)
    h.Settings.sep_odd_ind = {[1 2],[1 2],[1 2]};
    % for each set, ensure a number of leading standards 
    h.Settings.std_lead = [0 0 0]; % for each CP condition
    % min number of oddballs within each CP condition
    h.Settings.n_odd = [18, 18, 18]; % overrides h.Settings.totdur
    % number of minimum sets to randomise together
    h.Settings.n_set = [3 3 3]; % 1 = use min set size, for each CP condition
    % min number of oddballs per randomised set, per CP and per oddball
    % type
    h.Settings.n_odd_set = [6, 6, 6]; % overrides h.Settings.totdur
    % randomise sets?
    h.Settings.rand_set = 1; 
    
    %% RESPONSE PARAMETERS
    % record responses during experiment? 0 or 1
    h.Settings.record_response = 1;
    % how to record responses?
    h.Settings.record_response_type = {'all'}; %options: 'all','thistrial','previoustrial'
    % buttonpress options: key: keyboard inputs. Blank for no button press
    h.Settings.buttontype='key';
    % range of keyboard presses indicating a recordable response
    h.Settings.buttonopt = {'DownArrow'}; 
    
    %% THRESHOLDING
    % starting level and step size
    %h.Settings.threshold.startinglevel = 2; % for intensity)
    %h.Settings.threshold.step = 2;

    case 'NLT_classical'

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
    % names of any audiofiles
    h.Settings.audiofile = {};
    
    %% Condition-independent stimulus parameters - can be superceded by condition-dependent parameters
    % duration of stimulus sequence in seconds
    h.Settings.totdur = 0; 
    % duration of trial in seconds
    h.Settings.trialdur = 3; % if 0, consecutive stimuli will occur with no gap
    % duration of stimulus in seconds
    h.Settings.stimdur = 0; % modified by oddball settings
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
    h.Settings.inten = []; % value between 2 and 1000mA for Digitimer DS8R. leave blank if determined from GUI, threshold or adaptive
    % Tactile: number of pulses per trial
    h.Settings.nstim_trial = 1; % set to zero to be determined by stimdur
    % Tactile: within-trial frequency (Hz) 
    h.Settings.t_freq=0; 
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, intensityindex, pitch, channel
    h.Settings.oddballmethod = 'intensityindex'; % can use same type for pattern only if oddball intensity is adaptive
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
    % overridden by h.Settings.oddballvalue if using intensityindex
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
    h.Settings.n_set = [3 3 3]; % 1 = use min set size, for each CP condition
    % min number of oddballs within each CP condition
    h.Settings.n_odd = [36, 36, 36]; % overrides h.Settings.totdur
    % min number of oddballs per randomised set, per CP
    h.Settings.n_odd_set = [12, 12, 12]; % overrides h.Settings.totdur
    % randomise sets?
    h.Settings.rand_set = 1; 
    
    %% RESPONSE PARAMETERS
    % record responses during experiment? 0 or 1
    h.Settings.record_response = 1;
    % how to record responses?
    h.Settings.record_response_type = {'all'}; %options: 'all','thistrial','previoustrial'
    % buttonpress options: key: keyboard inputs. Blank for no button press
    h.Settings.buttontype='key';
    % range of keyboard presses indicating a recordable response
    h.Settings.buttonopt = {'DownArrow'}; 
    
    %% THRESHOLDING
    % starting level and step size
    %h.Settings.threshold.startinglevel = 2; % for intensity)
    %h.Settings.threshold.step = 2;
    
    case 'Complex'

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
    % names of any audiofiles
    h.Settings.audiofile = {};
    
    %% Condition-independent stimulus parameters - can be superceded by condition-dependent parameters
    % duration of stimulus sequence in seconds
    h.Settings.totdur = 0; 
    % duration of trial in seconds
    h.Settings.trialdur = 3; % if 0, consecutive stimuli will occur with no gap
    % duration of stimulus in seconds
    h.Settings.stimdur = 0; % modified by oddball settings
    % Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.patternmethod = 'intensity';
    h.Settings.patternvalue = [10 0 10 0 10 0 10 0 10 0]; % one per stimdur
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
    h.Settings.inten = []; % value between 2 and 1000mA for Digitimer DS8R
    % Tactile: number of pulses in a train
    
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
    h.Settings.ntrials = 50;
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
    h.Settings.buttonopt = {'DownArrow'}; 
    
    %% THRESHOLDING
    % starting level and step size
    %h.Settings.threshold.startinglevel = 2; % for intensity)
    %h.Settings.threshold.step = 2;
end

function h = setgeneral(h)

%% EQUIPMENT CONTROL
% record EEG, NS: netstation, BV: brainvision, 'serial': serial port
% serial port
h.Settings.serial = 'COM1';
%h.Settings.record_EEG='serial';
%h.Settings.EEGport = 'COM3'; % only needed for 'serial' EEG triggers
%h.Settings.EEGMarkPattern = 1; % mark EEG for every change in stimulus pattern (0 = start of trial only)
h.Settings.labjack=1; % Use labjack for controlling any equipment?
h.Settings.stimcontrol='LJTick-DAQ'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
h.Settings.labjack_DACport = 4;
h.Settings.stimchan = [6]; h.Settings.stimchanforLJ = 1;
h.Settings.nrchannels = 1; % total number of channels, e.g. on sound card
% channels on stimulator to use; use differenr rows for different pairs
% (e.g. for different conditions). If labjack, this can refer to output
% port(s).

% Tactile: number of pulses in a train
h.Settings.npulses_train = 1; % set to zero to be determined by stimdur
% Tactile: within-train frequency (Hz)
h.Settings.p_freq=[];
% Use timer to control frequency of labjack outputs? Otherwise uses software timing.
h.Settings.labjack_timer=1; % must use if there are a large number of pulses (e.g. >4)

% Audio sampling rate 
h.Settings.fs = 96000; % don't change this