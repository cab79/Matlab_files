function h = AdaptStair(h)

adaptive=0;
opt='';

if isfield(h,'s')
    s=h.s;
else
    s = AdaptStairParameters(h); % see function at the end
end

if ~isfield(s,'lasttrialrun')
    s.lasttrialrun = 0;
end

 % identify BUTTON PRESSES or OMISSIONS over previous trials of
% the same type
if h.Settings.adaptive.oddonly
    
    % only continue if next trial (that can be modified! i.e. h.i+h.Settings.ntrialsahead-1) is an oddball
    try
        if h.i+h.Settings.ntrialsahead-1 < length(h.Seq.signal)
            if any(ismember(h.Seq.odd_ind,h.Seq.condnum(h.i+h.Settings.ntrialsahead-1+1)))
                disp(['next (modifiable) trial is an oddball: trial ' num2str(h.i+h.Settings.ntrialsahead-1+1) ', condnum ' num2str(h.Seq.condnum(h.i+h.Settings.ntrialsahead-1+1))])
            end
            if h.Settings.adaptive.oddonly && ~any(ismember(h.Seq.odd_ind,h.Seq.condnum(h.i+h.Settings.ntrialsahead-1+1))) % h.odd_ind is generated in CreateSequence
                return
            end
        else
            return
        end
    catch
        error('define h.Seq.odd_ind and h.Seq.condnum in Sequence creation function');
    end

    % any responses over previous trials of this stimtype?
    notthis = find(h.Seq.signal~=h.Seq.signal(h.i)); % trials that are not this stimtype
    prevnotthis = notthis(notthis<h.i);% previous trials that were not this stimtype
    lastnotthis = max(prevnotthis); % last trial
    trls = lastnotthis+1:max(lastnotthis+h.Settings.adaptive.resptrials,h.i); % count to h.i or the max number of allowed trials
    trls_pressed = h.out.presstrial(ismember(h.out.presstrial,trls));
    % if there has been a button press that has not be
    % adapted to yet:
    if ~isempty(trls_pressed)
        h.pressedsinceoddball = 1;
    else
        h.omittedsinceoddball = 1;
    end
    
    % number of oddballs so far
    s.mintrialcount = sum(ismember(h.Seq.condnum(1:h.i),h.Seq.odd_ind));
else
    s.mintrialcount = h.i;
end

% run Adaptive if there have been BUTTON PRESSES on
% this trial or over previous trials that not yet been adapted to
if h.Settings.adaptive.oddonly
    if h.pressedsinceoddball && s.lasttrialrun~=trls_pressed(end)
        opt = 'responded';
        adaptive=1;
    end
end 
if (h.pressedlasttrial || h.pressedthistrial) && s.lasttrialrun~=h.i
    opt = 'responded';
    adaptive=1;
end

% if stimulation requires adaptive tuning to OMISSION of responses
if h.Settings.adaptive.omit && s.mintrialcount>h.Settings.adaptive.startomit
    if h.Settings.adaptive.oddonly
        if h.omittedsinceoddball
            opt = 'omitted';
            adaptive=1;
        end
    elseif h.omittedlasttrial && s.lasttrialrun~=h.i
        opt = 'omitted';
        adaptive=1;
    end
end

if ~adaptive
    return
end

% adapt to omitted response?
%omit=0;
%if ~isempty(varargin)
%   if strcmp(varargin{1},'omitted');
%       omit=1;
%   end
%end

%% RUN
disp('Running Adaptive')

% only run once per trial
%if strcmp(opt,'responded')
    if ~(s.lasttrialrun == h.i)
        s.lasttrialrun = h.i;
    else
        return
    end
%else
%    s.lasttrialrun = h.i;
%end

if size(s.expplan,1)==s.count_of_n_of_reversals
    s.expplan(end+1,:) = s.expplan(end,:);
end
    
% evaluate the subject's response
if strcmp(opt,'omitted')
    s.SubjectAccuracy(s.adaptive.trial)= 0;
elseif strcmp(opt,'responded')
    if h.pressedsinceoddball 
        presstrial=trls_pressed(end);
        pressbutton = h.out.pressbutton(h.out.presstrial==trls_pressed(end));
        correctsignal = h.Seq.signal(trls_pressed(end));
    elseif h.pressedlasttrial
        presstrial=h.i;
        pressbutton = h.out.lastpress{h.i};
        correctsignal = h.Seq.signal(h.i);
    elseif h.pressedthistrial
        presstrial=h.i;
        pressbutton = h.out.pressbutton(h.out.presstrial==h.i);
        correctsignal = h.Seq.signal(h.i);
    end
    
    resi = find(strcmp(pressbutton(end),h.Settings.buttonopt)); % which button was pressed?
    resfun = h.Settings.adaptive.signalval(resi); %what is meaning of this response?

    if resfun == correctsignal
    %s.SubjectAccuracy(s.adaptive.trial)= EvaluateAnswer(CorrectAnswer,s.feedback,Question);   %evaluing the subject answer (right or wrong)
        s.SubjectAccuracy(s.adaptive.trial)= 1;
    else
        s.SubjectAccuracy(s.adaptive.trial)= 0;
    end
end

% UPDATE THE ROWOFOUTPUT
s.rowofoutput (1, 1) = s.block;
s.rowofoutput (1, 2) = s.adaptive.trial;
s.rowofoutput (1, 3) = s.StimulusLevel;
s.rowofoutput (1, 4) = s.SubjectAccuracy(s.adaptive.trial);
% here I upddate the count for the up-s.down motion
if s.SubjectAccuracy(s.adaptive.trial) == 1
    s.n_down = s.n_down + 1;
    if s.n_down == s.down
        s.n_down = 0;
        s.pos = 1;
        s.trend = 1;
        % here I update the count of the number of s.reversals and
        % corresponding stepsize
        if s.pos ==1 && s.neg == -1
            s.count_of_n_of_reversals = s.count_of_n_of_reversals + 1;
            % calculate the threshold
            s.blockthresholds(s.n_threshold)=(s.StimulusLevel + s.rowofoutput(1, 3))/2;
            s.n_threshold = s.n_threshold + 1;
            s.actualstep = s.expplan(s.count_of_n_of_reversals, 2);
            s.pos = s.trend;
            s.neg = s.trend;
        end
        if s.isstep == 1
            s.StimulusLevel = s.StimulusLevel - s.actualstep;
        else
            s.StimulusLevel = s.StimulusLevel / s.actualstep;
        end
    end
else
    %error(['stopped at mintrialcount: ' num2str(s.mintrialcount)])
    s.neg = -1;
    s.trend = -1;
    s.n_down = 0;
    % here I update the count of the number of s.reversals and
    % corresponding stepsize
    if s.pos ==1 && s.neg == -1
        s.count_of_n_of_reversals = s.count_of_n_of_reversals + 1;
        % calculate the threshold
        s.blockthresholds(s.n_threshold)=(s.StimulusLevel + s.rowofoutput(1, 3))/2;
        s.n_threshold = s.n_threshold + 1;
        s.actualstep = s.expplan(s.count_of_n_of_reversals, 2);
        s.pos = s.trend;
        s.neg = s.trend;
    end
    if s.isstep == 1
        s.StimulusLevel = s.StimulusLevel + s.actualstep;
    else
        s.StimulusLevel = s.StimulusLevel * s.actualstep;
    end
end

% UPDATE THE ROWOFOUTPUT
s.rowofoutput (1, 5) = s.count_of_n_of_reversals;
s.rowofoutput (1, 6) = s.actualstep;
% update the number of trials
s.adaptive.trial = s.adaptive.trial + 1;

%disp(['length_blockthresh = ' num2str(length(s.blockthresholds))]);
disp(['nreversals = ' num2str(s.count_of_n_of_reversals)]);
disp(['level = ' num2str(s.StimulusLevel)]);

% here I calculate the threshol for the block
if length(s.blockthresholds)>=s.reversalForthresh
    switch s.thresholdtype
        case 'Arithmetic'
            s.expthresholds(s.block)=mean(s.blockthresholds(end-(s.reversalForthresh-1):end));
        case 'Geometric'
            s.expthresholds(s.block)=prod(s.blockthresholds(end-(s.reversalForthresh-1):end))^(1/length(s.blockthresholds(end-(s.reversalForthresh-1):end)));
        case 'Median'
            s.expthresholds(s.block)=median(s.blockthresholds(end-(s.reversalForthresh-1):end));
        otherwise
            disp('Unknown calculation type.')
    end
    fprintf('Threshold equal to %1.3f\n', s.expthresholds(s.block));
else
    s.expthresholds(s.block)=NaN;
end
s.rowofoutput (1, 7) = s.expthresholds(s.block);

% UPDATE THE GLOBAL OUTPUT VARIABLE
s.out.adaptive = [s.out.adaptive; s.rowofoutput];
h.s =s;
h.out.adaptive = s.out.adaptive;
%fprintf('Press return to continue\n');
%pause
%fprintf ('\nBLOCK ENDED\n');
%pause(2)

function s = AdaptStairParameters(h)

% REMOVE THESE PARAMETERS LATER

%s.STIMULUSLEVEL = INIZIALIZES A VECTOR FOR THE LEVEL OF THE STIMULI
%s.SUBJCTACCURACY = INIZIALIZES A VECTOR FOR EVALUATION OF SUBJECT ACCURACY
%s.FA = INIZIALIZES A VECTOR FOR THE FALSE ALLARMS
%s.FEEDBACK= LOGICAL VALUE. DISPLAY ANSWER CORRECTNESS
%s.TEMPORARYTHRESHOLD = TEMPORARY THRESHOLDS AFTER EACH TRIAL
%s.MATSAVEDATA = INIZIALIZES A MATRIX FOR SAVING THE DATA

%s.tracking = 'Staircase';
s.block=1;
%s.standard= -30; % attentuation applied to the standard tone - needs to be less that max possible intensity
s.feature='TwoDownOneUp';
s.down=2;
s.feedback= 1;
s.reversals = [4;8];
s.stepsize = [2;sqrt(2)];
%s.fileout = 'data.txt';
s.SaveResults =1;
s.tasktype = 1;
s.exppos = 1;
%s.experiment = 'IntensityDiscriminationPureTone_S';
s.thresholdtype='Arithmetic';
s.reversalForthresh = 3;
s.isstep = 0;
s.nafc=3;
%s.NameStepSize='Factor';

% if this is the first run, do some setup
if ~isfield(s,'count_of_n_of_reversals')
    if length(s.reversals) ~= length(s.stepsize)
        error('The number of s.reversals and the number of steps must be identical.');
    end
    % here I set the plan of the threshold tracking, i.e., a matrix (two
    % columns) that contains on the left the progressive number of s.reversals and
    % on the right the corresponding step size
    s.expplan = [(1:sum(s.reversals))', zeros(sum(s.reversals), 1)];
    i=1;
    for j=1:length(s.reversals)
        for k=1:s.reversals(j)
            s.expplan(i, 2)=s.stepsize(j);
            i=i+1;
        end
    end
    % here I define the variable 'row of output' that contains all the output values of the
    % function. In the current function the output is updated at the end of the while loop
    % this are the values and the labels
    s.out.adaptive = [];
    s.rowofoutput = zeros(1, 6);
    s.expthresholds = zeros(1, 1);

    %clc
    %input('Press return to begin ', 's');
    %pause(1)
    % indexes for the while loop
    s.count_of_n_of_reversals = 0;
    s.adaptive.trial = 1;
    s.blockthresholds = zeros(length(s.reversalForthresh), 1);
    s.n_threshold = 1;
    % variable for the up-s.down
    s.n_down = 0;
    % variable for count the positive and negative answers
    s.pos = 0;
    s.neg = 0;
    s.trend = 30;
    s.StimulusLevel = h.Settings.adaptive.startinglevel;
    s.actualstep = s.expplan(1, 2);
end