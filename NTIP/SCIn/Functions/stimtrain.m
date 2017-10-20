function h = stimtrain(h,opt)

% mods needed: use settings to define port number


switch h.Settings.stimcontrol

    case 'labjack'
        
        if ~isfield(h,'ljHandle')
            try
                h.ljHandle = get(h.ljhandle, 'Value');
            end
        end
        
        ljud_LoadDriver; % Loads LabJack UD Function Library
        ljud_Constants; % Loads LabJack UD constant file
        if h.Settings.stimchanforLJ
            port = h.Settings.stimchan;
        else
            port=4;
        end
        
        switch opt
            case 'start'

                if h.Settings.labjack_timer
                    %Set the timer/counter pin offset to 4, which will put the first timer/counter on FIO4. 
                    Error = ljud_AddRequest (h.ljHandle,  LJ_ioPUT_CONFIG, LJ_chTIMER_COUNTER_PIN_OFFSET, port, 0, 0);
                    Error_Message(Error)

                    %use 48MHz clock base with divisor = 48 to get 1 MHz timer clock: 
                    Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_CONFIG, LJ_chTIMER_CLOCK_BASE, LJ_tc1MHZ_DIV, 0, 0);
                    Error_Message(Error)
                    Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_CONFIG, LJ_chTIMER_CLOCK_DIVISOR, 200, 0, 0); 
                    Error_Message(Error)

                    %Enable 2 timers.
                    Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_CONFIG, LJ_chNUMBER_TIMERS_ENABLED, 2, 0, 0); 
                    Error_Message(Error)

                    %Configure Timer0 as Frequency out.
                    Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_TIMER_MODE, 0, LJ_tmFREQOUT, 0, 0); 
                    Error_Message(Error)

                    %Set the second divisor to N (x2), yielding a frequency of 1000000/(N*2) Hz
                    N = 1000000 / (h.Settings.p_freq*2);
                    Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_TIMER_VALUE, 0, 250, 0, 0); 
                    Error_Message(Error)
                    
                    %Configure Timer1 as timer stop:
                    Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_TIMER_MODE, 1, LJ_tmTIMERSTOP, 0, 0);
                    Error_Message(Error)

                    %set number of pulses: 
                    Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_TIMER_VALUE, 1, 10, 0, 0);
                    Error_Message(Error)

                    %Execute the requests. 
                    Error = ljud_GoOne(h.ljHandle);
                    Error_Message(Error)
                    disp('running')
                else
                    % pulse train instruction
                    for pr = 1:h.Settings.npulses_train % train
                        Error = ljud_AddRequest(h.ljHandle,LJ_ioPUT_DIGITAL_BIT,port,1,0,0); % 
                        Error_Message(Error)

                        Error = ljud_AddRequest(h.ljHandle,LJ_ioPUT_WAIT,port,round((1000000/h.Settings.p_freq)/2),0,0); % Actual resolution is 64 microseconds.
                        Error_Message(Error)

                        Error = ljud_AddRequest(h.ljHandle,LJ_ioPUT_DIGITAL_BIT,port,0,0,0);
                        Error_Message(Error)

                        Error = ljud_AddRequest(h.ljHandle,LJ_ioPUT_WAIT,port,round((1000000/h.Settings.p_freq)/2),0,0); % Actual resolution is 64 microseconds.
                        Error_Message(Error)
                    end
                    %Execute the stimulus train
                    t1=GetSecs;
                    Error = ljud_GoOne(h.ljHandle);
                    Error_Message(Error)
                    %ljud_GetResult(ljHandle, LJ_ioGET_DIGITAL_BIT, 7, @Value)
                    t2=GetSecs;
                    disp(['pr: ' num2str(pr) '; train ISI: ' num2str(h.Settings.p_freq) 's;  labjack stimulus: ' num2str(t2-t1)]);
                end
                
            case 'stop'
                
                if h.Settings.labjack_timer
                    Error = ljud_AddRequest (h.ljHandle,  LJ_ioPUT_CONFIG, LJ_chTIMER_COUNTER_PIN_OFFSET, port, 0, 0);
                    Error_Message(Error)
                    %Execute the requests. 
                    Error = ljud_GoOne(h.ljHandle);
                    Error_Message(Error)
                end
         end
        
    case 'audioplayer'
        if ~exist('opt','var')
            opt = 'run';
        end
        
        switch opt
            case 'run'
                h=sinwave(h);
                h.Seq.aud = audioplayer(h.Seq.stimseq', h.Settings.fs);
                play(h.Seq.aud);
                %sound(h.Seq.stimseq', h.Settings.fs, 16);
                %pause(0.2)
            
            case 'create'
                h=sinwave(h);
                h.Seq.aud = audioplayer(h.Seq.stimseq', h.Settings.fs);
            
            case 'start'
                play(h.Seq.aud);
                h.playstart = GetSecs;
                
            case 'pause'
                pause(h.Seq.aud);

            case 'resume'
                resume(h.Seq.aud);

            case 'stop'
                stop(h.Seq.aud);
                
            case 'getsample'
                h.currentsample=get(h.Seq.aud,'CurrentSample');
                h.totalsamples=get(h.Seq.aud,'TotalSamples');
        end
        
    case 'PsychPortAudio'
        switch opt
            case 'setup'
                h = PTBaudio(h);
                
            case 'getsample'
                s = PsychPortAudio('GetStatus', h.pahandle);
                if s.Active == 1
                    h.currentsample=s.ElapsedOutSamples;
                    %h.totalsamples=;
                end
                
            case 'create' 
                h=sinwave(h);
                %if strcmp(h.Settings.design,'trials')
                %    PsychPortAudio('FillBuffer', h.pahandle, h.Seq.stimseq);
                %    
                %elseif strcmp(h.Settings.design,'continuous')
                %    h.pabuffer = PsychPortAudio('CreateBuffer', h.pahandle, h.Seq.stimseq);% Engine still running on a schedule?
                %   
                %end
                
            case 'start' 
                if strcmp(h.Settings.design,'trials')
                    PsychPortAudio('FillBuffer', h.pahandle, h.Seq.stimseq);
                    PsychPortAudio('Start', h.pahandle, 1, 0, 1);

                elseif strcmp(h.Settings.design,'continuous')
                    h.pabuffer = PsychPortAudio('CreateBuffer', h.pahandle, h.Seq.stimseq);
                   
                    s = PsychPortAudio('GetStatus', h.pahandle);
                    if s.Active == 0 %&& ~isfield(h,'i') % new run
                        PsychPortAudio('UseSchedule', h.pahandle, 1, length(h.Seq.signal));
                        PsychPortAudio('AddToSchedule', h.pahandle, h.pabuffer);
                        h.playstart = PsychPortAudio('Start', h.pahandle, 0, 0, 1);
                    %elseif s.Active == 0 
                    %    error('increase ntrialsahead in Settings')
                    else
                        PsychPortAudio('AddToSchedule', h.pahandle, h.pabuffer);
                        disp('new trial(s) added to schedule')
                    end
                end
        end
end
