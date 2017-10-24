function h = Experiment(h,opt)

% GUI handle name
h.GUIhname = findall(0, 'Type', 'figure', 'Tag', 'SCIn');

switch opt
    
    case 'setup'
        
        if isfield(h.Settings,'buttontype')
            if ~isempty(h.Settings.buttontype)
                opt = 'setup';
                h = buttonpress(h,opt);
            end
        end
        
        if isfield(h.Settings,'labjack')
            if h.Settings.labjack
                opt = 'connect';
                h = Labjack(h,opt);
            end
        end
        
        if isfield(h.Settings,'record_EEG')
            if h.Settings.record_EEG
                opt = 'connect';
                h = recordEEG(h,opt);
            end
        end
        
        if isfield(h.Settings,'D188')
            if h.Settings.D188
                opt = 'setup';
                h = D188(h,opt);
            end
        end
        
        %if isfield(h.Settings,'stimcontrol') 
        %    if strcmp(h.Settings.stimcontrol,'PsychPortAudio')
        %        h = PTBaudio(h);
        %    end
        %end
        
        % setup stim control
        if isfield(h.Settings,'stimcontrol')
            if ~isempty(h.Settings.stimcontrol)
                opt = 'setup';
                h = stimtrain(h,opt);
            end
        end
        
        % create all trials if design is continuous
        if isfield(h.Settings,'stimcontrol') && strcmp(h.Settings.design,'continuous') && ~isfield(h.Seq,'stimseq')
            if ~isempty(h.Settings.stimcontrol)
                opt = 'create';
                h = stimtrain(h,opt);
            end
        end
    
    case 'start'
        
        % create output structure
        allout = cell(1,length(h.Seq.signal));
        h.out.stimtime = allout;
        
        % note time that experiment started
        h.t_start = datestr(now,30);
        
        % record stim sequence
        if isfield(h,'Seq')
            global d
            [~,seqname,~] = fileparts(h.SeqName);
            fname = ['Sequence_' h.subID '_' seqname '_startblock' num2str(h.startblock) '_' h.t_start '.mat'];
            seq = h.Seq;
            save(fullfile(d.root,d.out,fname),'seq');
        end
          
        % KBCheck or KBQueueCheck results
        KbName('UnifyKeyNames'); %used for cross-platform compatibility of keynaming
        h.out.presstrial = [];
        h.out.pressbutton = [];
        h.out.presstime = [];
        h.out.RT = [];
        h.out.presstimedelta = [];

        % KBQueueCheck results
        h.out.firstpressbutton = allout;
        h.out.firstpress = allout;
        h.out.firstrelease = allout;
        h.out.lastpress = allout;
        h.out.lastrelease = allout;
        
        if isfield(h.Settings,'buttontype')
            if ~isempty(h.Settings.buttontype)
                if h.Settings.record_response
                    KbQueueCreate; %creates cue using defaults
                    KbQueueStart;  %starts the cue
                end
            end
        end
        
        % require input to start, e.g. scanner trigger, msgbox
        if isfield(h.Settings,'blockstart')
            if ~isempty(h.Settings.blockstart)
                if iscell(h.Settings.blockstart)
                    for i = 1:length(h.Settings.blockstart)
                        opt = h.Settings.blockstart{i};
                        h = blockstart(h,opt);
                    end
                else
                    opt = h.Settings.blockstart;
                    h = blockstart(h,opt);
                end
            end
        end
        
        if isfield(h.Settings,'record_EEG')
            if h.Settings.record_EEG
                opt = 'start';
                h = recordEEG(h,opt);
            end
        end
        
        % record scanner triggers
        %if isfield(h,'trigtime')
        %    h.out.trigtime = h.trigtime;
        %end
        
        
        % Flush Buffer so only responses after stim onset are recorded (not recording accidental buton presses)
        KbQueueFlush; 
        
        switch h.Settings.design
            case 'trials'
                h = SeqLoop(h);
            case 'continuous'
                h = SeqRun(h);
        end
        
    case {'pause','resume','stop'}
        if isfield(h.Settings,'stimcontrol')
            if ~isempty(h.Settings.stimcontrol)
                h = stimtrain(h,opt); % stimulus train
            end
        end
    
end

function h = SeqLoop(h)
Priority(2);
tic

h.i=0;

while h.i<length(h.Seq.signal)
    
    % record previous start time to calculate ISI
    if isfield(h,'st')
        h.st_prev = h.st; 
    end
    
    % start time of trial
    h.st = GetSecs;
    h.ct=h.st;
    
    % calculate ISI
    if isfield(h,'st_prev')
        isi = h.st-h.st_prev;
    else
        isi = 0;
    end
    
    t=toc/60;
    h.i=h.i+1;
    disp(['Trial' num2str(h.i) '. Elapsed time is ' num2str(t) ' mins. ISI is ' num2str(isi) ' s']);
    
    % D188 - set output channel
    if isfield(h,'D188')
        %if h.Settings.D188
            opt = 'setchan';
            h.D188.chan = h.Seq.signal(h.i);
            h = D188(h,opt);
        %end
    end
    
    % send stimulus
    if isfield(h.Settings,'stimcontrol')
        if ~isempty(h.Settings.stimcontrol)
            opt = 'create';
            h = stimtrain(h,opt); % stimulus train
            opt = 'start';
            h = stimtrain(h,opt); % stimulus train
        end
    end
    
    % record stimulus timing
    h.out.stimtime{h.i} = GetSecs;
    
    % STIM marker on EEG
    if isfield(h.Settings,'record_EEG')
        if h.Settings.record_EEG
            opt = 'mark';
            h = recordEEG(h,opt);
        end
    end
    
    
    % define the duration of this trial
    h.trialdur = h.Settings.trialdur; % if stimdur is zero, will start next trial immediately after stimulus ends
    
    % enter waiting time loop 
    h.stop=0;
    h.savedRT=0;
    h.savedi=0;
    h = waitingloop(h);
    
    if h.stop==1
        break
    end
    
    % record first and last responses on previous trial
    h = record_response(h,'prev_trial');
end

% for continuous sequences
function h = SeqRun(h)
Priority(2);
tic

% send stimulus
if isfield(h.Settings,'stimcontrol')
    if ~isempty(h.Settings.stimcontrol)
        if h.Settings.ntrialsahead
            opt = 'create';
            h = stimtrain(h,opt); % stimulus train
        end
        opt = 'start';
        h = stimtrain(h,opt); % stimulus train
        h.laststimload=GetSecs;
        opt = 'getsample';
        h = stimtrain(h,opt);
    end
end

% trial duration
if h.Settings.ntrialsahead
    if isfield(h,'totalsamples')
        h.trialdur = h.totalsamples/h.Settings.fs;
    else
        try
            h.trialdur = h.totdur;
        catch
            h.trialdur = sum(h.Settings.stimdur)*length(h.Seq.signal);
        end
    end
else
    h.trialdur = size(h.Seq.stimseq,2);
end

% start time of trial
h.st = GetSecs;
h.ct=h.st;

% wait and record responses 
h.stop=0;
h.savedRT=0;
h.savedi=0;
h.i=1;
h.out.stimtime{h.i} = h.st;
h.marktime = cell(1,length(h.Seq.signal));
disp(['Trial ' num2str(h.i) '/' num2str(length(h.Seq.signal)) '.']);
h = waitingloop(h);

function h = waitingloop(h)

global d
h.ct=GetSecs;
while (h.ct-h.st)<h.trialdur
    
    %WaitSecs(0.001);
        %h.ct=GetSecs;disp(['exit time = ' num2str(h.ct-h.st)]);
        
    
    % update current time and exit if needed
    if strcmp(h.Settings.design,'trials')
        h.ct=GetSecs;
        if (h.ct-h.st)>h.trialdur
            break; 
        end
    end
    
%t2=GetSecs;
%try
%    disp(['delay = ' num2str(t2-t1)]);
%end

    % run special function for continuous sequences
    if strcmp(h.Settings.design,'continuous')
        h = ContFun(h);
    end
    h = record_response(h,'all');

    %only continue loop (which can be slow) if not too close to end of
    %trial. Increases accuracy of EEG markers.
    t1=GetSecs;
    if h.i>1 && h.i<length(h.Seq.signal)
        if h.out.projend{h.i}>t1+0.1
            continue
        end
    end
    
    %
    
    % get GUI data for pause/resume
    GUIh = guihandles(h.GUIhname);

    % update current time and exit if needed
    if strcmp(h.Settings.design,'trials')
        h.ct=GetSecs;
        if (h.ct-h.st)>h.trialdur; 
            break; 
        end;
    end
    
    try
        px = get(GUIh.PauseResume, 'Value');
    catch
        px=0;
    end
    if ~isfield(h,'disableGUI')
        pause(0.001);
    end
    while px
        
        % pause stimulus
        if isfield(h.Settings,'stimcontrol')
            if ~isempty(h.Settings.stimcontrol)
                opt = 'pause';
                h = stimtrain(h,opt);
            end
        end
        
        % Mark EEG
        if isfield(h.Settings,'record_EEG')
            if h.Settings.record_EEG
                opt = 'pause';
                h = recordEEG(h,opt);
            end
        end

        % get GUI data for pause/resume
        %GUIh = guihandles(h.GUIhname);
        try
            rx = ~get(GUIh.PauseResume, 'Value');
        catch
            rx=0;
        end
        if ~isfield(h,'disableGUI')
            pause(0.001);
        end
        if rx   
            
            % Mark EEG
            if isfield(h.Settings,'record_EEG')
                if h.Settings.record_EEG
                    opt = 'resume';
                    h = recordEEG(h,opt);
                end
            end

            % option to re-start on button press
            % require input to start, e.g. scanner trigger, msgbox
            if isfield(h.Settings,'blockstart')
                if ~isempty(h.Settings.blockstart)
                    if iscell(h.Settings.blockstart)
                        for i = 1:length(h.Settings.blockstart)
                            opt = h.Settings.blockstart{i};
                            h = blockstart(h,opt);
                        end
                    else
                        opt = h.Settings.blockstart;
                        h = blockstart(h,opt);
                    end
                end
            end
            
            % resume stimulus
            if isfield(h.Settings,'stimcontrol')
                if ~isempty(h.Settings.stimcontrol)
                    opt = 'resume';
                    h = stimtrain(h,opt);
                end
            end
            
            
            % update trial duration
            if h.Settings.ntrialsahead
                if isfield(h,'totalsamples')
                    h.trialdur = (h.totalsamples-h.currentsample)/h.Settings.fs;
                else
                    try
                        h.trialdur = h.totdur - (GetSecs-h.st);
                    catch
                        h.trialdur = sum(h.Settings.stimdur)*length(h.Seq.signal) - (GetSecs-h.st);
                    end
                end
            else
                h.trialdur = size(h.Seq.stimseq,2) - h.i;
            end
            
            % start time of trial
            h.st = GetSecs;
            h.ct=h.st;

            px=0;
        end
    end

    % update current time and exit if needed
    if strcmp(h.Settings.design,'trials')
        h.ct=GetSecs;
        if (h.ct-h.st)>h.trialdur; 
            break; 
        end;
    end

    % get GUI data for start/stop
    %GUIh = guihandles(h.GUIhname);
    try
        sx = ~get(GUIh.StartStop, 'Value');
    catch
        sx=0;
    end
    if ~isfield(h,'disableGUI')
        pause(0.001);
    end
    while sx
        
        % stop stimulus
        if isfield(h.Settings,'stimcontrol')
            if ~isempty(h.Settings.stimcontrol)
                opt = 'stop';
                h = stimtrain(h,opt);
            end
        end
        
        h.stop = 1;
        if isfield(h,'pahandle')
            PsychPortAudio('Close', h.pahandle);
        end
        
        try
            global spt1
            fclose(spt1);
        end
       
        
        % release keyboard cue
        KbQueueRelease; 
        
        break
        %pause(0.1)
        %while get(StartStop, 'Value')        
        %    break;
        %end
    end
    if h.stop==1
        break
    end

    % update current time and exit if needed
    if strcmp(h.Settings.design,'trials')
        h.ct=GetSecs;
        if (h.ct-h.st)>h.trialdur; 
            break; 
        end;
    end

    if h.i~=length(h.Seq.signal)
        % next trial is in a new block, pause here
        if h.Seq.blocks(h.i+1)>h.Seq.blocks(h.i)
            
            % Mark EEG
            if isfield(h.Settings,'record_EEG')
                if h.Settings.record_EEG
                    opt = 'pause';
                    h = recordEEG(h,opt);
                end
            end
            
            % require input to start, e.g. scanner trigger, msgbox
            if isfield(h.Settings,'blockstart')
                if ~isempty(h.Settings.blockstart)
                    if iscell(h.Settings.blockstart)
                        for i = 1:length(h.Settings.blockstart)
                            opt = h.Settings.blockstart{i};
                            h = blockstart(h,opt);
                        end
                    else
                        opt = h.Settings.blockstart;
                        h = blockstart(h,opt);
                    end
                end
            end
            
            % Mark EEG
            if isfield(h.Settings,'record_EEG')
                if h.Settings.record_EEG
                    opt = 'resume';
                    h = recordEEG(h,opt);
                end
            end
            tic
        end
    else
        disp('END OF EXPERIMENT');
        break
    end

    % update current time and exit if needed
    if strcmp(h.Settings.design,'trials')
        h.ct=GetSecs;
        if (h.ct-h.st)>h.trialdur; 
            break; 
        end;
    end

    %if ~isempty(h.out.RT{h.i}) %&& h.Settings.record_response==1 %&& h.trialdur-(h.ct-h.st) > 0.1 % only save if 0.1s left, otherwise maybe too slow
    if h.savedi ~=h.i % save on each new trial, but no later
        [~,seqname,~] = fileparts(h.SeqName);
        fname = ['Output_' h.subID '_' seqname '_startblock' num2str(h.startblock) '_' h.t_start '.mat'];
        out = h.out;
        save(fullfile(d.root,d.out,fname),'out');

        % monitor what has just been saved so as not to save it again
        % (speeds up the loop)
        %h.savedRT =h.out.RT{h.i};
        h.savedi =h.i;
    end
    
    % update current time and exit if needed
    if strcmp(h.Settings.design,'trials')
        h.ct=GetSecs;
        if (h.ct-h.st)>h.trialdur
            break; 
        end
    end
end

function h = ContFun(h)
% find current sample of stim
if isfield(h.Settings,'stimcontrol')
    if ~isempty(h.Settings.stimcontrol)
        opt = 'getsample';
        h = stimtrain(h,opt);
    end
end
nowtime = GetSecs;
%disp([num2str(h.currentsample) ' / ' num2str(h.totalsamples)])

% mark EEG when pattern changes occur mid-trial
if isfield(h.Settings,'record_EEG')
    if h.Settings.record_EEG
        if h.Settings.EEGMarkPattern
            
            if isempty(h.marktime{h.i})
                % on this trial, what time to mark?
                h.marktime{h.i} = repmat(h.out.stimtime{h.i},1,length(h.Seq.PatternSamples{h.i,1})) + h.Seq.PatternSamples{h.i,1}/h.Settings.fs;

                % don't mark the end of the last event as this is a new trial (marked later)
                if ~isempty(h.marktime{h.i})
                    h.marktime{h.i}(end) = [];
                end

                % number of within-trial marks needed on this trial
                h.imark{h.i} = 1:length(h.marktime{h.i});
            end

            % has the first time been exceeded?
            if ~isempty(h.imark{h.i})
                if h.marktime{h.i}(h.imark{h.i}(1))<nowtime
                    % STIM marker on EEG
                    opt = 'mark';
                    h = recordEEG(h,opt);

                    % remove marker index just marked
                    if ~isempty(h.imark{h.i})
                        h.imark{h.i}(1) = [];
                    end
                end
            end
        end
    end
end

% projected end time of trial
if h.i>1
    triallength = h.Seq.trialend(h.i)-h.Seq.trialend(h.i-1);
    proj_end = h.out.stimtime{h.i} + triallength/h.Settings.fs;
elseif h.i==1
    triallength = h.Seq.trialend(h.i);
    proj_end = h.st + triallength/h.Settings.fs;
%else
%    triallength = 0;
%    proj_end = h.st + h.Seq.trialend(1)/h.Settings.fs;
end

% new-trial events:
% if current sample is greater than the sample at the end of the trial 
% OR if projected time of the end sample is greater than the current time
trials = find(h.Seq.trialend > h.currentsample);
newtrial=0;
if ~isempty(trials)
    if h.i ~= trials(1)
        newtrial = 1;
    end
end
if nowtime>=proj_end
    newtrial = 1;
end

if newtrial
    %try
    %    disp(['Currentsample error: ' num2str((h.currentsample - h.Seq.trialend(h.i))/h.Settings.fs)])
    %catch
    %    disp(['Currentsample error: ' num2str(h.currentsample - h.Seq.trialend(1))])
    %end
    %disp(['Current time error: ' num2str(nowtime - proj_end)])
    
    % record first and last responses on previous trial
    h = record_response(h,'prev_trial');
    
    %update h.i
    %try
        %h.i = h.i+1;
    try
        h.i = max(h.i+1,trials(1));
    catch
        h.i = h.i +1;
    end
    %catch
    %    h.i = 1;
    %end
    
    % record stimulus timing
    h.out.stimtime{h.i} = nowtime;
    
    % ISI
    try
        h.out.isi{h.i-1} = h.out.stimtime{h.i}-h.out.stimtime{h.i-1};

        % expected isi
        h.out.expisi{h.i-1} = triallength/h.Settings.fs;
        % discrepency
        h.out.discrep{h.i-1} = h.out.isi{h.i-1} - h.out.expisi{h.i-1};
        
        % projected end of this trial in s (for speeding up loop)
        
        triallength = (h.Seq.trialend(h.i)-h.Seq.trialend(h.i-1))/h.Settings.fs;
        h.out.projend{h.i} = h.out.stimtime{h.i}+triallength;
    end
    
    % D188 - set output channel
    if isfield(h,'D188')
        %if h.Settings.D188
            opt = 'setchan';
            h.D188.chan = h.Seq.signal(h.i);
            h = D188(h,opt);
        %end
    end
    
    % STIM marker on EEG
    if isfield(h.Settings,'record_EEG')
        if h.Settings.record_EEG
            opt = 'mark';
            h = recordEEG(h,opt);
        end
    end
    
    % add next n trials to the buffer (unless buffer full)
    % send stimulus
    if isfield(h.Settings,'stimcontrol') && h.i<=length(h.Seq.signal)-h.Settings.ntrialsahead && h.Settings.ntrialsahead
        if ~isempty(h.Settings.stimcontrol)
            opt = 'create';
            h = stimtrain(h,opt); % stimulus train
            opt = 'start';
            h = stimtrain(h,opt); % stimulus train
            h.laststimload = GetSecs;
        end
    end
    
    % display
    %t=toc/60;
    try
    %    tot = h.totalsamples/h.Settings.fs/60;
    %    disp(['Trial ' num2str(h.i) '. Elapsed time is ' num2str(t) '/' num2str(tot) ' mins. ISI is ' num2str(isi) ' s. Onset discrepency: ' num2str(h.out.discrep{h.i})]);
    %catch
        disp(['Trial ' num2str(h.i) '/' num2str(length(h.Seq.signal)) '. Channel ' num2str(h.Seq.signal(h.i)) '. Last ISI was ' num2str(h.out.isi{h.i-1}) ' s. ']);
    end
    
    %if num2str(h.out.isi{h.i-1})<0.9
    %    dbstop in Experiment.m at 686
    %end
    
    if h.Settings.plottrials
        h = plot_wave(h,h.Settings.plottrials);
    end
    
end

function h = record_response(h,opt)

% initialise
h.pressedlasttrial = 0;
h.omittedlasttrial = 0;
h.pressedthistrial = 0;
h.pressedsinceoddball = 0;
h.omittedsinceoddball = 0;

switch opt
   case 'prev_trial'
        % record first and last responses on this trial
        if isfield(h.Settings,'buttontype')
            if ~isempty(h.Settings.buttontype)
                if h.Settings.record_response
                    [h.pressedlasttrial, firstPress, firstRelease, lastPress, lastRelease] = KbQueueCheck;
                    if h.pressedlasttrial
                        h.out.firstpressbutton{h.i} = KbName(firstPress);
                        h.out.firstpress{h.i} = firstPress(firstPress>0);
                        h.out.firstrelease{h.i} = firstRelease(firstRelease>0);
                        h.out.lastpressbutton{h.i} = KbName(lastPress);
                        h.out.lastpress{h.i} = lastPress(lastPress>0);
                        h.out.lastrelease{h.i} = lastRelease(lastRelease>0);
                    else
                        h.omittedlasttrial=1;
                    end
                end
            end
        end
   case 'current_trial'
        if isfield(h.Settings,'buttontype')
            if ~isempty(h.Settings.buttontype)
                if h.Settings.record_response

                    % OPTION 1: get GUI data for button press
                  %  GUIh = guihandles(h.GUIhname);
                  %  if isfield(GUIh,'buttonpressed')
                  %      btnstr=get(GUIh.buttonpressed,'String');
                  %      presstime = get(GUIh.buttontime,'String');
                  %      if ~strcmp(btnstr,'Waiting')

                    % OPTION 2: get kbcheck data
                    [keyIsDown,presstime, keyCode, deltaSecs] = KbCheck;
                    btnstr = KbName(keyCode);

                    if keyIsDown && ~isempty(btnstr)
                        h.pressedthistrial=1;
                        recordresp=1;
                        % if button options are specified
                        if isfield(h.Settings,'buttonopt')
                            if ~isempty(h.Settings.buttonopt)
                                if ~ismember(btnstr,h.Settings.buttonopt)
                                    recordresp=0;
                                end
                            end
                        end
                        if recordresp
                            h.out.presstrial = [h.out.presstrial h.i];
                            h.out.pressbutton = [h.out.pressbutton {btnstr}];
                            h.out.presstime = [h.out.presstime presstime];
                            h.out.presstimedelta = [h.out.presstimedelta deltaSecs];
                            if strcmp(h.Settings.design,'trials')
                                h.out.RT = [h.out.RT presstime-h.st];
                            elseif strcmp(h.Settings.design,'continuous')
                                h.out.RT = [h.out.RT presstime-h.out.stimtime{h.i}];
                            end
                            disp(['RESPONSE: ' btnstr])
                        end


                        % clear button press buffer if in continuous mode
                        % (trial mode clears at the beginning of each new trial)
                        %if isfield(h.Settings,'buttontype') %&& strcmp(h.Settings.design,'continuous')
                        %    if ~isempty(h.Settings.buttontype)
                        %        opt = 'setup';
                        %        h = buttonpress(h,opt);
                        %    end
                        %end
                    end
                end
            end
        end
        
     case 'all'
        if isfield(h.Settings,'buttontype')
            if ~isempty(h.Settings.buttontype)
                if h.Settings.record_response
                    [h.pressed, firstPress, firstRelease, lastPress, lastRelease] = KbQueueCheck;
                    if h.pressed
                        recordresp=1;
                        fpress = KbName(firstPress);
                        lpress = KbName(lastPress);
                        firstPress = firstPress(firstPress>0);
                        lastPress = lastPress(lastPress>0);
                        % if button options are specified
                        if isfield(h.Settings,'buttonopt')
                            if ~isempty(h.Settings.buttonopt)
                                if ~ismember(fpress,h.Settings.buttonopt) || ~ismember(lpress,h.Settings.buttonopt)
                                    recordresp=0;
                                end
                            end
                        end
                        if recordresp
                            if firstPress~=lastPress
                                h.out.presstrial = [h.out.presstrial h.i*ones(1, length(firstPress)) h.i*ones(1, length(lastPress))];
                                h.out.pressbutton = [h.out.pressbutton fpress lpress];
                                h.out.presstime = [h.out.presstime firstPress lastPress];
                                h.out.presstimedelta = [h.out.presstimedelta zeros(1, length(firstPress)) zeros(1, length(lastPress))];
                                if strcmp(h.Settings.design,'trials')
                                    h.out.RT = [h.out.RT firstPress-h.st lastPress-h.st];
                                elseif strcmp(h.Settings.design,'continuous')
                                    h.out.RT = [h.out.RT firstPress-h.out.stimtime{h.i} lastPress-h.out.stimtime{h.i}];
                                end
                                disp(['RESPONSE: ' fpress ' ' lpress])
                            else
                                h.out.presstrial = [h.out.presstrial h.i*ones(1, length(firstPress))];
                                h.out.pressbutton = [h.out.pressbutton fpress];
                                h.out.presstime = [h.out.presstime firstPress];
                                h.out.presstimedelta = [h.out.presstimedelta zeros(1, length(firstPress))];
                                if strcmp(h.Settings.design,'trials')
                                    h.out.RT = [h.out.RT firstPress-h.st];
                                elseif strcmp(h.Settings.design,'continuous')
                                    h.out.RT = [h.out.RT firstPress-h.out.stimtime{h.i}];
                                end
                                disp(['RESPONSE: ' fpress])
                            end
                            % Flush Buffer so only new responses are recorded next
                            % time
                            KbQueueFlush; 
                        end
                           
                    end
                    
                end
            end
        end
end

% if stimulation requires adaptive tuning to responses
if isfield(h.Settings,'adaptive')
    if ~isempty(h.Settings.adaptive)
        h = AdaptStair(h);
    end
end

function h = plot_wave(h,ntrials)

if ~isfield(h,'f')
    h.f = figure; 
end

try
    samples = h.Seq.trialend(h.i-ntrials):h.Seq.trialend(h.i);
    %samples=samples-round(length(samples)/ntrials);
    sig = h.Seq.stimseq(:,samples(samples>1));

    figure(h.f)
    subplot(3,1,1)
    plot(1/96000:1/96000:size(sig,2)/96000,sig(1,:))
    subplot(3,1,2)
    plot(1/96000:1/96000:size(sig,2)/96000,sig(2,:))
    subplot(3,1,3)
    plot(1/96000:1/96000:size(sig,2)/96000,sig(1,:)-sig(2,:))
    title(['Trial ' num2str(h.i)])
    %if h.i==21
    %    pause
    %end
end