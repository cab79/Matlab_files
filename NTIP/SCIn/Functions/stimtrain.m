function h = stimtrain(h,opt)

switch h.Settings.stim(h.trialstimnum).control

    case {'labjack','LJTick-DAQ'}
        
        switch opt
            case 'create'
                h=trial(h,'set&construct');
                return
        end
        
        if ~isfield(h,'ljHandle')
            try
                h.ljHandle = get(h.ljhandle, 'Value');
            end
        end
        
        ljud_LoadDriver; % Loads LabJack UD Function Library
        ljud_Constants; % Loads LabJack UD constant file
        if h.Settings.stim(h.trialstimnum).chanforLJ
            port = h.Settings.stim(h.trialstimnum).chan;
        else
            port=6;
        end
        
        switch opt
            case 'setup'
                % if GUI is set to intensity > 0, check this is not
                % accidental
                if isfield(h,'inten_mean_gui'); 
                    inten_mean = str2double(h.inten_mean_gui);
                    if inten_mean>0
                        choice = questdlg(['Is ' num2str(inten_mean) ' the correct starting intensity?'], ...
                        'Intensity', ...
                        'Yes','No','Yes');
                        % Handle response
                        switch choice
                            case 'No'
                                error('wrong intensity set in GUI')
                        end
                    end
                end
            
            case 'create'
                h=trial(h,'set&construct');
                
            case 'calc'
                h=trial(h,'set');
                
            case 'getsample'
                h.currentsample = round((h.ct-h.out.stimtime{1})*h.Settings.fs);
            
            case 'set'
                
                % if selecting the intensity from a time sample of a
                % waveform
                if strcmp(h.Settings.design,'continuous')
                    if isfield(h,'currentsample')
                        if h.currentsample>0
                            h.inten = h.Seq.stimseq(h.currentsample);
                            h.out.stimseq_record(h.currentsample) = h.inten;
                        else
                            h.inten = 0;
                        end
                    else
                        h.inten = 0;
                    end
                end
                
                
                %Set DACA 
                if strcmp(h.Settings.stim(h.trialstimnum).control,'LJTick-DAQ')
                    try
                        Error = ljud_ePut(h.ljHandle, LJ_ioTDAC_COMMUNICATION, LJ_chTDAC_UPDATE_DACA, h.inten*h.Settings.DAC_multiply, 0); 
                        Error_Message(Error)
                    catch% to use DAC0 port
                        Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_DAC, 0, h.inten*h.Settings.DAC_multiply, 0,0);
                        Error_Message(Error)
                        Error = ljud_GoOne(h.ljHandle);
                        Error_Message(Error)
                    end
                elseif strcmp(h.Settings.stim(h.trialstimnum).control,'labjack')
                    Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_DAC, 0, h.inten*h.Settings.DAC_multiply, 0,0);
                    Error_Message(Error)
                    Error = ljud_GoOne(h.ljHandle);
                    Error_Message(Error)
                end

                %voltage = 80;
                %current = h.inten*100;
                %time = 50/1e6 * h.Settings.npulses_train;
                %mJ = voltage * current * time
               
                
            case 'start'

                if h.Settings.stim(h.trialstimnum).labjack_timer && ~isempty(h.Settings.stim(h.trialstimnum).p_freq)
                    if h.Settings.stim(h.trialstimnum).p_freq>=h.LJfreqtable(1,1)
                        %Set the timer/counter pin offset to 6, which will put the first timer/counter on FIO6. 
                        Error = ljud_AddRequest (h.ljHandle,  LJ_ioPUT_CONFIG, LJ_chTIMER_COUNTER_PIN_OFFSET, port, 0, 0);
                        Error_Message(Error)

                        % get row of table. Columns are Hz, base clock, clock divisor, and timer value. 
                        r = dsearchn(h.LJfreqtable(:,1),h.Settings.stim(h.trialstimnum).p_freq);
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
                        div = round(12e6 / (h.Settings.stim(h.trialstimnum).p_freq * 2^16));
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
                    Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_TIMER_VALUE, 1, h.Settings.stim(h.trialstimnum).npulses_train, 0, 0);
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
                    %port=2
                    if isempty(h.Settings.stim(h.trialstimnum).p_freq) || h.Settings.stim(h.trialstimnum).p_freq == 0
                       p_freq = 10000; % 0.1 ms delay (DS8R can detect down to 0.01ms)
                    else
                       p_freq = h.Settings.stim(h.trialstimnum).p_freq; 
                    end
                    for pr = 1:h.Settings.stim(h.trialstimnum).npulses_train % train
                        Error = ljud_AddRequest(h.ljHandle,LJ_ioPUT_DIGITAL_BIT,port,1,0,0); % 
                        Error_Message(Error)

                        Error = ljud_AddRequest(h.ljHandle,LJ_ioPUT_WAIT,port,round((1000000/p_freq)/2),0,0); % Actual resolution is 64 microseconds.
                        Error_Message(Error)

                        Error = ljud_AddRequest(h.ljHandle,LJ_ioPUT_DIGITAL_BIT,port,0,0,0);
                        Error_Message(Error)

                        Error = ljud_AddRequest(h.ljHandle,LJ_ioPUT_WAIT,port,round((1000000/p_freq)/2),0,0); % Actual resolution is 64 microseconds.
                        Error_Message(Error)
                    end
                    %Execute the stimulus train
                    t1=GetSecs;
                    Error = ljud_GoOne(h.ljHandle);
                    Error_Message(Error)
                    %ljud_GetResult(ljHandle, LJ_ioGET_DIGITAL_BIT, 7, @Value)
                    t2=GetSecs;
                    disp(['pulses per stim: ' num2str(pr) ';  labjack stim length: ' num2str(t2-t1)]);
                end
                
            case 'stop'
                
                if strcmp(h.Settings.stim(h.trialstimnum).control,'LJTick-DAQ')
                    Error = ljud_ePut(h.ljHandle, LJ_ioTDAC_COMMUNICATION, LJ_chTDAC_UPDATE_DACA, 0, 0); 
                    Error_Message(Error)
                else
                    %try
                        Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_DAC, 0, 0, 0,0);
                        Error_Message(Error)
                        Error = ljud_GoOne(h.ljHandle);
                        Error_Message(Error)
                    %end
                end
                
                if isfield(h.Settings.stim(h.trialstimnum),'labjack_timer')
                    if h.Settings.stim(h.trialstimnum).labjack_timer
                        Error = ljud_AddRequest (h.ljHandle,  LJ_ioPUT_CONFIG, LJ_chTIMER_COUNTER_PIN_OFFSET, port, 0, 0);
                        Error_Message(Error)
                        %Execute the requests. 
                        Error = ljud_GoOne(h.ljHandle);
                        Error_Message(Error)
                    end
                end
                
        end
         
    case 'serial'
        switch opt
            case 'setup'
                opt = 'spt1';
                open_serial(h,opt);
                opt = 'spt';
                open_serial(h,opt);
            case 'set'
                h=trial(h,'set');
                global spt
                if (h.inten~=0 && h.inten<=255)
                    fprintf(spt,'%s', h.inten);
                else
                    error('invalid Intensity level')
                end

            case 'start'
                % trigger stimulator and mark EEG at the same time
                global spt1
                TriggerNum = h.Seq.signal(h.trialstimnum,h.i);
                fprintf(spt1,num2str(32+TriggerNum));
        end

    case 'audioplayer'
        if ~exist('opt','var')
            opt = 'run';
        end
        
        switch opt
            case 'run'
                h=trial(h,'set&construct');
                h.Seq.aud = audioplayer(h.Seq.stimseq', h.Settings.fs);
                play(h.Seq.aud);
                %sound(h.Seq.stimseq', h.Settings.fs, 16);
                %pause(0.2)
            
            case 'create'
                h=trial(h,'set&construct');
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
                h=trial(h,'set&construct');
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
                        PsychPortAudio('UseSchedule', h.pahandle, 1, size(h.Seq.signal,2));
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

