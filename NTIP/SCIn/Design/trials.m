function h = trials(h,opt)

% GUI handle name
h.GUIhname = findall(0, 'Type', 'figure', 'Tag', 'SCIn');

%get(h.GUIhname, 'WindowKeyPressFcn')
%set(h.GUIhname, 'WindowKeyPressFcn', []);
%get(h.GUIhname, 'WindowKeyPressFcn')
%set(h.GUIhname, 'WindowKeyPressFcn', @(hObject,eventdata) SCIn('SCIn_WindowKeyPressFcn',hObject,eventdata,guidata(hObject)));
%get(h.GUIhname, 'WindowKeyPressFcn')
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
        
        if isfield(h.Settings,'stimcontrol') 
            if strcmp(h.Settings.stimcontrol,'PsychPortAudio')
                opt = 'setup';
                h = PTBaudio(h,opt);
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
        
        if isfield(h.Settings,'blockstart')
            if ~isempty(h.Settings.blockstart)
                h = blockstart(h);
            end
        end
        
        if isfield(h.Settings,'buttontype')
            if ~isempty(h.Settings.buttontype)
                if h.Settings.buttonstart
                    opt = 'start';
                    h = buttonpress(h,opt);
                end
            end
        end
        
        if isfield(h.Settings,'record_EEG')
            if h.Settings.record_EEG
                opt = 'start';
                h = recordEEG(h,opt);
            end
        end
        
        % create output structure
        allout = cell(1,length(h.Seq.signal));
        h.out.stimtime = allout;
        
        if isfield(h.Settings,'buttontype')
            if ~isempty(h.Settings.buttontype)
                if h.Settings.record_response
                    KbName('UnifyKeyNames'); %used for cross-platform compatibility of keynaming
                    KbQueueCreate; %creates cue using defaults
                    KbQueueStart;  %starts the cue
                    
                    % KBCheck results
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
                end
            end
        end
        
        % note time that experiment started
        h.t_start = datestr(now,30);
        
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

while h.i<length(h.Seq.signal);
    
    % record previous start time to calculate ISI
    if isfield(h,'st'); 
        h.st_prev = h.st; 
    end;
    
    % start time of trial
    h.st = GetSecs;
    h.ct=h.st;
    
    % calculate ISI
    if isfield(h,'st_prev');
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
            h = stimtrain(h); % stimulus train
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
    
    % Flush Buffer so only responses after stim onset are recorded
    KbQueueFlush; 
    
    % define the duration of this trial
    h.trialdur = 1/h.Settings.freq;
    
    %set(h.GUIhname, 'WindowKeyPressFcn', @(hObject,eventdata) SCIn('SCIn_WindowKeyPressFcn',hObject,eventdata,guidata(hObject)));
                %pause(0.001);
               
    
    % enter waiting time loop 
    h.stop=0;
    h.savedRT=0;
    h.savedi=0;
    h = waitingloop(h);
    
    if h.stop==1
        break
    end
    
    % record first and last responses on this trial
    if isfield(h.Settings,'buttontype')
        if ~isempty(h.Settings.buttontype)
            if h.Settings.record_response
                [pressed, firstPress, firstRelease, lastPress, lastRelease] = KbQueueCheck;
                if pressed
                    h.out.firstpressbutton{h.i} = KbName(firstPress);
                    h.out.firstpress{h.i} = firstPress(firstPress>0);
                    h.out.firstrelease{h.i} = firstRelease(firstRelease>0);
                    h.out.lastpress{h.i} = lastPress(lastPress>0);
                    h.out.lastrelease{h.i} = lastRelease(lastRelease>0);
                end
            end
        end
    end
end

% for continuous sequences
function h = SeqRun(h)
Priority(2);
tic

% send stimulus
if isfield(h.Settings,'stimcontrol')
    if ~isempty(h.Settings.stimcontrol)
        opt = 'start';
        h = stimtrain(h,opt); % stimulus train
        opt = 'getsample';
        h = stimtrain(h,opt); % stimulus train
    end
end

% trial duration
h.trialdur = h.totalsamples/h.Settings.fs;

% start time of trial
h.st = GetSecs;
h.ct=h.st;

% wait and record responses 
h.stop=0;
h.savedRT=0;
h.savedi=0;
h.i=0;
h = waitingloop(h);

function h = waitingloop(h)

global d
h.ct=GetSecs;
while (h.ct-h.st)<h.trialdur
    
    pause(0.001);
        %h.ct=GetSecs;disp(['exit time = ' num2str(h.ct-h.st)]);
    
    % update current time and exit if needed
    if strcmp(h.Settings.design,'trials')
        h.ct=GetSecs;
        if (h.ct-h.st)>h.trialdur; 
            break; 
        end;
    end

    % run special function for continuous sequences
    if strcmp(h.Settings.design,'continuous')
        h = ContFun(h);
    end
    
%t2=GetSecs;
%try
%    disp(['delay = ' num2str(t2-t1)]);
%end
%t1=GetSecs; 

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
                
                    if keyIsDown
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
                                h.out.RT = [h.out.RT presstime-h.out.expstimtime{h.i}];
                            end
                            disp(['RESPONSE: ' btnstr])
                        end
                        
                        % if stimulation requires adaptive tuning to responses
                        if isfield(h.Settings,'adaptive')
                            if ~isempty(h.Settings.adaptive)
                                h = AdaptStair(h);
                            end
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
                %end
            end
        end
    end
    
    
    % update current time and exit if needed
    if strcmp(h.Settings.design,'trials')
        h.ct=GetSecs;
        if (h.ct-h.st)>h.trialdur; 
            break; 
        end;
    end
    
    % get GUI data for pause/resume
    GUIh = guihandles(h.GUIhname);
    try
        px = get(GUIh.PauseResume, 'Value');
    catch
        px=0;
    end
    pause(0.001);
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
        GUIh = guihandles(h.GUIhname);
        try
            rx = ~get(GUIh.PauseResume, 'Value');
        catch
            rx=0;
        end
        pause(0.001);
        if rx   
            
            % Mark EEG
            if isfield(h.Settings,'record_EEG')
                if h.Settings.record_EEG
                    opt = 'resume';
                    h = recordEEG(h,opt);
                end
            end

            % option to re-start on button press
            if isfield(h.Settings,'buttontype')
                if ~isempty(h.Settings.buttontype)
                    if h.Settings.buttonstart
                        opt = 'start';
                        h = buttonpress(h,opt);
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
    GUIh = guihandles(h.GUIhname);
    try
        sx = ~get(GUIh.StartStop, 'Value');
    catch
        sx=0;
    end
    pause(0.001);
    while sx
        
        % stop stimulus
        if isfield(h.Settings,'stimcontrol')
            if ~isempty(h.Settings.stimcontrol)
                opt = 'pause';
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
            
            h = blockstart(h);
            
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
        if (h.ct-h.st)>h.trialdur; 
            break; 
        end;
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
%disp([num2str(h.currentsample) ' / ' num2str(h.totalsamples)])

% trial events
trials = find(h.Seq.trialend > h.currentsample);
if h.i ~= trials(1)
    
    % record first and last responses on previous trial
    if h.i>1
        if isfield(h.Settings,'buttontype')
            if ~isempty(h.Settings.buttontype)
                if h.Settings.record_response
                    [pressed, firstPress, firstRelease, lastPress, lastRelease] = KbQueueCheck;
                    if pressed
                        h.out.firstpressbutton{h.i} = KbName(firstPress);
                        h.out.firstpress{h.i} = firstPress(firstPress>0);
                        h.out.firstrelease{h.i} = firstRelease(firstRelease>0);
                        h.out.lastpress{h.i} = lastPress(lastPress>0);
                        h.out.lastrelease{h.i} = lastRelease(lastRelease>0);
                    end
                end
            end
        end
    end
    
    %create h.i
    h.i = trials(1);
    
    % ISI
    try
        isi = (h.Seq.trialend(h.i) - h.Seq.trialend(h.i-1))/h.Settings.fs;
    catch
        isi = h.Seq.trialend(h.i)/h.Settings.fs;
    end
    
    % record stimulus timing
    h.out.stimtime{h.i} = GetSecs;
    % expected timing
    h.out.expstimtime{h.i} = h.playstart + (h.Seq.trialend(h.i)/h.Settings.fs - isi);
    % discrepency
    h.out.discrep{h.i} = h.out.expstimtime{h.i} - h.out.stimtime{h.i};
    
    % STIM marker on EEG
    if isfield(h.Settings,'record_EEG')
        if h.Settings.record_EEG
            opt = 'mark';
            h = recordEEG(h,opt);
        end
    end
    
    % Flush Buffer so only responses after stim onset are recorded
    KbQueueFlush; 
    
    % display
    t=toc/60;
    tot = h.totalsamples/h.Settings.fs/60;
    disp(['Trial' num2str(h.i) '. Elapsed time is ' num2str(t) '/' num2str(tot) ' mins. ISI is ' num2str(isi) ' s. Onset discrepency: ' num2str(h.out.discrep{h.i})]);
    % send EEG trigger on new trial
    
end
