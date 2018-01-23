function h = stimtrain(h,opt)


switch h.Settings.stimcontrol

    case {'labjack','LJTick-DAQ'}
        
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
            port=6;
        end
        
        switch opt
            case 'setDAC'
                
                %convert h.Settings.inten from mA to V
                inten = h.Settings.inten/100;
                
                if length(inten)==1
                    % if THRESHOLD
                    if isfield(h.Settings,'threshold')
                        if ~isempty(h.Settings.threshold)
                            if ~isfield(h,'s')
                                inten = h.Settings.threshold.startinglevel/100;
                            else
                                inten = h.s.actStimulusLevel/100;
                            end
                        end
                    end

                    %Set DACA 
                    Error = ljud_ePut(h.ljHandle, LJ_ioTDAC_COMMUNICATION, LJ_chTDAC_UPDATE_DACA, inten, 0); 
                    Error_Message(Error)
                elseif isfield(h.Settings,'p_freq') % apply a train of intensity changes at this freq
                    for i = 1:length(inten)
                        %Set DACA 
                        Error = ljud_AddRequest(h.ljHandle, LJ_ioTDAC_COMMUNICATION, LJ_chTDAC_UPDATE_DACA, inten, 0, 0); 
                        Error_Message(Error)
                    
                        Error = ljud_AddRequest(h.ljHandle,LJ_ioPUT_WAIT,h.Settings.labjack_DACport,round(1000000/h.Settings.p_freq),0,0); % Actual resolution is 64 microseconds.
                        Error_Message(Error)
                    end
                end

                % to use DAC0 port
                %Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_DAC, 0, 0.00, 0,0);
                %Error_Message(Error)
                
                Error = ljud_GoOne(h.ljHandle);
                Error_Message(Error)
               
                
            case 'start'
                
                % projected time at end of trial
                h.out.projend{h.i} = h.ct+h.Settings.trialdur;

                if h.Settings.labjack_timer
                    if h.Settings.p_freq>=h.LJfreqtable(1,1)
                        %Set the timer/counter pin offset to 6, which will put the first timer/counter on FIO6. 
                        Error = ljud_AddRequest (h.ljHandle,  LJ_ioPUT_CONFIG, LJ_chTIMER_COUNTER_PIN_OFFSET, port, 0, 0);
                        Error_Message(Error)

                        % get row of table. Columns are Hz, base clock, clock divisor, and timer value. 
                        r = dsearchn(h.LJfreqtable(:,1),h.Settings.p_freq);
                        % get parameters
                        freq = h.LJfreqtable(r,1);
                        base = h.LJfreqtable(r,2);
                        div = h.LJfreqtable(r,3);
                        val = h.LJfreqtable(r,4);
                        disp(['actual freq is ' num2str(freq)]);
                        
                        if base == 1e6; basenum = 23;%   //1 MHz clock base w/ divisor (no Counter0)
                        elseif base == 4e6; basenum = 24;%   //4 MHz clock base w/ divisor (no Counter0)
                        elseif base == 12e6; basenum = 25;%  //12 MHz clock base w/ divisor (no Counter0)
                        elseif base == 48e6; basenum = 26;%  //48 MHz clock base w/ divisor (no Counter0)
                        end

                        %use 48MHz clock base with divisor = 48 to get 1 MHz timer clock: 
                        Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_CONFIG, LJ_chTIMER_CLOCK_BASE, basenum, 0, 0);
                        Error_Message(Error)
                        Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_CONFIG, LJ_chTIMER_CLOCK_DIVISOR, div, 0, 0); 
                        Error_Message(Error)

                        %Enable 2 timers.
                        Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_CONFIG, LJ_chNUMBER_TIMERS_ENABLED, 2, 0, 0); 
                        Error_Message(Error)

                        %Configure Timer0 as Frequency out.
                        Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_TIMER_MODE, 0, LJ_tmFREQOUT, 0, 0); 
                        Error_Message(Error)

                        Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_TIMER_VALUE, 0, val, 0, 0); 
                        Error_Message(Error)
                        
                    else % use less accurate method for low freq
                %Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_DAC, 0, 0.5, 0,0);
                %Error_Message(Error)
                %Error = ljud_GoOne(h.ljHandle);
                %Error_Message(Error)
                    
                        %Set the timer/counter pin offset to 6, which will put the first timer/counter on FIO6. 
                        Error = ljud_AddRequest (h.ljHandle,  LJ_ioPUT_CONFIG, LJ_chTIMER_COUNTER_PIN_OFFSET, port, 0, 0);
                        Error_Message(Error)

                        %Enable 2 timers.
                        Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_CONFIG, LJ_chNUMBER_TIMERS_ENABLED, 2, 0, 0); 
                        Error_Message(Error)

                        % use 12MHz clock base
                        Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_CONFIG, LJ_chTIMER_CLOCK_BASE, LJ_tc12MHZ_DIV, 0, 0);
                        Error_Message(Error)

                        % control freq output with a divisor from 6 to 183 to get an output frequency of about 30.5Hz to 1Hz respectively: 
                        div = round(12e6 / (h.Settings.p_freq * 2^16));
                        freq = 12e6/div/2^16;
                        disp(['actual freq is ' num2str(freq)]);
                        %div = round(12e6 / (0.1 * 2^16));
                        Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_CONFIG, LJ_chTIMER_CLOCK_DIVISOR, div, 0, 0); 
                        Error_Message(Error)

                        %Configure Timer0 as 16-bit PWM.
                        Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_TIMER_MODE, 0, LJ_tmPWM16, 0, 0); 
                        Error_Message(Error)

                        %Initialize the 16-bit PWM with a 50% duty cycle.
                        Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_TIMER_VALUE, 0, 32767, 0, 0); 
                        Error_Message(Error)
                    end

                    %Configure Timer1 as timer stop:
                    Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_TIMER_MODE, 1, LJ_tmTIMERSTOP, 0, 0);
                    Error_Message(Error)

                    %set number of pulses: 
                    Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_TIMER_VALUE, 1, h.Settings.npulses_train, 0, 0);
                    Error_Message(Error)

                    %Execute the requests. 
                    Error = ljud_GoOne(h.ljHandle);
                    Error_Message(Error)
                    disp('running')
                    
                    
                %Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_DAC, 0, 0, 0,0);
                %Error_Message(Error)
                %Error = ljud_GoOne(h.ljHandle);
                %Error_Message(Error)
                    
                else
                    % pulse train instruction
                    port=2
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
                
                try
                    Error = ljud_ePut(h.ljHandle, LJ_ioTDAC_COMMUNICATION, LJ_chTDAC_UPDATE_DACA, 0, 0); 
                    Error_Message(Error)
                end
                
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
