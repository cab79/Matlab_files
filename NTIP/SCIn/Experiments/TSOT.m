function h = TSOT(h,opt)

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
    
    case 'start'
        
        mb_handle = msgbox('Put electrodes on participant and press OK to continue.','Message');
        uiwait(mb_handle);
        
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
                    h.out.pressbutton = allout;
                    h.out.presstime = allout;
                    h.out.RT = allout;
                end
            end
        end
        
        % note time that experiment started
        h.t_start = datestr(now,30);
        
        %h.pause = 0;
        %h.stop = 0;
        h = SeqLoop(h);
        
    %case 'pause'
    %    h.pause=1;
    
    %case 'resume'
    %    h.pause=0;
    
    %case 'stop'
    %    stop(h.Seq.signal);
end

function h = SeqLoop(h)
Priority(2);
tic

h.i=0;
while h.i<=length(h.Seq.signal);
    
    % start time of trial
    h.st = GetSecs;
    h.ct=h.st;
    
    t=toc/60;
    h.i=h.i+1;
    disp(['Trial' num2str(h.i) '. Elapsed time is ' num2str(t) ' mins']);
    
    % D188 - set output channel
    if isfield(h.Settings,'D188')
        if h.Settings.D188
            opt = 'setchan';
            h.D188.chan = h.Seq(h.i);
            h = D188(h,opt);
        end
    end
    
    % send stimulus
    if isfield(h.Settings,'stimcontrol')
        if ~isempty(h.Settings.stimcontrol)
            h = stimtrain(h); % stimulus train
        end
    end
    
    % record stimulus timing
    h.out.stimtime(h.i) = GetSecs;
    
    % STIM marker on EEG
    if isfield(h.Settings,'record_EEG')
        if h.Settings.record_EEG
            opt = 'mark';
            h = recordEEG(h,opt);
        end
    end
    
    % define the duration of this trial
    h.trialdur = 1/h.Settings.freq;
    
    % enter waiting time loop 
    h = waitingloop(h);
    
end

function h = waitingloop(h)

global d
h.ct=GetSecs;
while h.ct-h.st>h.trialdur

    if isfield(h.Settings,'buttontype')
        if ~isempty(h.Settings.buttontype)
            if h.Settings.record_response
                % get GUI data for button press
                GUIh = guidata(h.GUIhname);
                if isfield(GUIh.button,'press')
                    recordresp=1;
                    % if button options are specified
                    if isfield(h.Settings,'buttonopt')
                        if ~isempty(h.Settings.buttonopt)
                            if ~ismember(GUIh.button.press,h.Settings.buttonopt)
                                recordresp=0;
                            end
                        end
                    end
                    if recordresp
                        h.out.pressbutton{h.i} = GUIh.button.press;
                        h.out.presstime{h.i} = h.button.presstime;
                        h.out.RT{h.i} = h.button.presstime-st;
                    end
                    % if stimulation requires adaptive tuning to responses
                    if isfield(h.Settings,'adaptive')
                        if ~isempty(h.Settings.adaptive)
                            
                        end
                    end
                    GUIh = rmfield(GUIh.button,'press');
                    guidata(h.GUIhname,GUIh)
                end
            end
        end
    end

    % update current time and exit if needed
    h.ct=GetSecs;
    if h.ct-h.st>h.trialdur; break; end;

    while get(GUIh.PauseResume, 'Value')
        pause(0.1)

        % Mark EEG
        if isfield(h.Settings,'record_EEG')
            if h.Settings.record_EEG
                opt = 'pause';
                h = recordEEG(h,opt);
            end
        end

        while ~get(GUIh.PauseResume, 'Value')   
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

            break;
        end
    end

    % update current time and exit if needed
    h.ct=GetSecs;
    if h.ct-h.st>h.trialdur; break; end;

    stop=0;
    while ~get(GUIh.StartStop, 'Value')
        stop = 1;
        break
        %pause(0.1)
        %while get(StartStop, 'Value')        
        %    break;
        %end
    end
    if stop==1
        break
    end

    % update current time and exit if needed
    h.ct=GetSecs;
    if h.ct-h.st>h.trialdur; break; end;

    if h.i~=length(h.Seq.signal)
        % next trial is in a new block, pause here
        if h.Seq.blocks(h.i+1)>h.Seq.blocks(h.i)
            mb_handle = msgbox('End of block: Press OK to start the next block.','Message');

            % Mark EEG
            if isfield(h.Settings,'record_EEG')
                if h.Settings.record_EEG
                    opt = 'pause';
                    h = recordEEG(h,opt);
                end
            end

            uiwait(mb_handle);
            pause(1);
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
    end

    % update current time and exit if needed
    h.ct=GetSecs;
    if h.ct-h.st>h.trialdur; break; end;

    if h.Settings.record_response==1
        fname = ['Output_' h.subID '_' h.SeqName '_startblock' num2str(h.startblock) '_' h.t_start];
        out = h.out;
        save(fullfile(d.root,d.out,fname),'out');
    end
end

