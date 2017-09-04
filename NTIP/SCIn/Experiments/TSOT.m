function h = TSOT(h,opt)

% GUI handle name
h.GUIhname = findall(0, 'Type', 'figure', 'Tag', 'SCIn');

switch opt
    
    case 'setup'
        
        BPopt = 'setup';
        h = buttonpress(h,BPopt);
        
        NSopt = 'connect';
        h = Netstation(h,NSopt);
        
        LJopt = 'connect';
        h = Labjack(h,LJopt);
        
        D1opt = 'setup';
        h = D188(h,D1opt);
    
    case 'start'
        
        mb_handle = msgbox('Put electrodes on participant and press OK to continue.','Message');
        uiwait(mb_handle);
        
        %if isfield(h,'pedal')
        %    FPopt = 'start';
        %    h = footpedal(h,FPopt);
        %end
        
        if isfield(h,'button')
            BPopt = 'start';
            h = buttonpress(h,BPopt);
        end
        
        if isfield(h,'NS')
            NSopt = 'start';
            h = Netstation(h,NSopt);
        end
        
        % create output structure
        allout = cell(1,length(h.Seq.signal));
        h.out.stimtime = allout;
        h.out.pressbutton = allout;
        h.out.presstime = allout;
        h.out.RT = allout;
        
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
global d
Priority(2);
tic

h.i=0;
while h.i<=length(h.Seq.signal);
    
    % start time of trial
    st = GetSecs;
    ct=st;
    
    t=toc/60;
    h.i=h.i+1;
    disp(['Trial' num2str(h.i) '. Elapsed time is ' num2str(t) ' mins']);
    
    % D188 - set output channel
    if isfield(h,'D188')
        D1opt = 'setchan';
        h.D188.chan = h.Seq(h.i);
        h = D188(h,D1opt);
    end
    
    % send stimulus
    h = stimtrain(h);
    
    % record stimulus timing
    h.out.stimtime(h.i) = GetSecs;
    
    % STIM marker on EEG
    if isfield(h,'NS') % Netstation
        NSopt = 'mark';
        h = Netstation(h,NSopt);
    elseif isfield(h,'BV') % Brainvision
        BVopt = 'mark';
        h = Brainvision(h,BVopt);
    end
    
    % update current time and exit if needed
    ct=GetSecs;
    if ct-st>h.freq; break; end;
    
    % get GUI data for button press
    GUIh = guidata(h.GUIhname);
    if isfield(GUIh.button,'press')
        %if ismember(GUIh.button.press,h.button.opt)
            h.out.pressbutton{h.i} = GUIh.button.press;
            h.out.presstime{h.i} = h.button.presstime;
            h.out.RT{h.i} = h.button.presstime-st;
        %end
        GUIh = rmfield(GUIh.button,'press');
        guidata(h.GUIhname,GUIh)
    else
        h.out.pressbutton{h.i} = NaN;
        h.out.presstime{h.i} = NaN;
        h.out.RT{h.i} = NaN;
    end
    
    % update current time and exit if needed
    ct=GetSecs;
    if ct-st>h.freq; break; end;
    
    while get(GUIh.PauseResume, 'Value')
        pause(0.1)
        
        % Mark EEG
        if isfield(h,'NS') % Netstation
            NSopt = 'pause';
            h = Netstation(h,NSopt);
        elseif isfield(h,'BV') % Brainvision
            BVopt = 'pause';
            h = Brainvision(h,BVopt);
        end
        
        while ~get(GUIh.PauseResume, 'Value')   
            % Mark EEG
            if isfield(h,'NS') % Netstation
                NSopt = 'resume';
                h = Netstation(h,NSopt);
            elseif isfield(h,'BV') % Brainvision
                BVopt = 'resume';
                h = Brainvision(h,BVopt);
            end
            
            % option to re-start on button press
            if isfield(h,'button')
                BPopt = 'start';
                h = buttonpress(h,BPopt);
            end
            
            break;
        end
    end
    
    % update current time and exit if needed
    ct=GetSecs;
    if ct-st>h.freq; break; end;
    
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
    ct=GetSecs;
    if ct-st>h.freq; break; end;
    
    if h.i~=length(h.Seq.signal)
        if h.Seq.blocks(h.i+1)>h.Seq.blocks(h.i)
            mb_handle = msgbox('End of block: Press OK to start the next block.','Message');
            
            % Mark EEG
            if isfield(h,'NS') % Netstation
                NSopt = 'pause';
                h = Netstation(h,NSopt);
            elseif isfield(h,'BV') % Brainvision
                BVopt = 'pause';
                h = Brainvision(h,BVopt);
            end
            
            uiwait(mb_handle);
            pause(1);
            % Mark EEG
            if isfield(h,'NS') % Netstation
                NSopt = 'resume';
                h = Netstation(h,NSopt);
            elseif isfield(h,'BV') % Brainvision
                BVopt = 'resume';
                h = Brainvision(h,BVopt);
            end  
            tic
        end
    else
        disp('END OF EXPERIMENT');
    end
    
    % update current time and exit if needed
    ct=GetSecs;
    if ct-st>h.freq; break; end;
    
    if h.record_response==1
        fname = ['Output_' h.subID '_' h.SeqName '_startblock' num2str(h.startblock) '_' t_start];
        out = h.out;
        save(fullfile(d.root,d.out,fname),'out');
    end
end

function h = stimtrain(h)

train_dur = round((1000/h.p_freq)*h.npulses_train);
% pulse train instruction
for pr = 1:h.npulses_train % train
    Error = ljud_AddRequest(ljHandle,LJ_ioPUT_DIGITAL_BIT,4,1,0,0);
    Error_Message(Error)

    Error = ljud_AddRequest(ljHandle,LJ_ioPUT_WAIT,4,round((1000000/h.p_freq)/2),0,0);
    Error_Message(Error)

    Error = ljud_AddRequest(ljHandle,LJ_ioPUT_DIGITAL_BIT,4,0,0,0);
    Error_Message(Error)

    Error = ljud_AddRequest(ljHandle,LJ_ioPUT_WAIT,4,round((1000000/h.p_freq)/2),0,0);
    Error_Message(Error)
end
%Execute the stimulus train
Error = ljud_GoOne(ljHandle);
Error_Message(Error)
%ljud_GetResult(ljHandle, LJ_ioGET_DIGITAL_BIT, 7, @Value)


