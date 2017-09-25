function h = stimtrain(h,opt)

% mods needed: use settings to define port number


switch h.Settings.stimcontrol

    case 'labjack'
        % pulse train instruction
        for pr = 1:h.Settings.npulses_train % train
            Error = ljud_AddRequest(ljHandle,LJ_ioPUT_DIGITAL_BIT,4,1,0,0);
            Error_Message(Error)

            Error = ljud_AddRequest(ljHandle,LJ_ioPUT_WAIT,4,round((1000000/h.Settings.p_freq)/2),0,0);
            Error_Message(Error)

            Error = ljud_AddRequest(ljHandle,LJ_ioPUT_DIGITAL_BIT,4,0,0,0);
            Error_Message(Error)

            Error = ljud_AddRequest(ljHandle,LJ_ioPUT_WAIT,4,round((1000000/h.Settings.p_freq)/2),0,0);
            Error_Message(Error)
        end
        %Execute the stimulus train
        Error = ljud_GoOne(ljHandle);
        Error_Message(Error)
        %ljud_GetResult(ljHandle, LJ_ioGET_DIGITAL_BIT, 7, @Value)
        
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
        h=sinwave(h);
        
        % Fill the audio playback buffer with the audio data:
        %[h.actualStartTime, ~, ~, h.estStopTime] = 
        PsychPortAudio('Stop', h.pahandle);%, 0, 1);
        PsychPortAudio('FillBuffer', h.pahandle, h.Seq.stimseq);
        %PsychPortAudio('Volume', h.pahandle,h.inten(1));
        
        % Compute new start time for follow-up beep, beepPauseTime after end of
        % previous one
        startCue = 0;%h.estStopTime + h.Settings.dur;

        % Start audio playback
        % Should we wait for the device to really start (1 = yes)
        % INFO: See help PsychPortAudio
        waitForDeviceStart = 0;
        PsychPortAudio('Start', h.pahandle, 1, startCue, waitForDeviceStart);
end
