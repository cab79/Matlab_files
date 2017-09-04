function h = Brainvision(h,opt)

switch opt
    case 'connect'
        %% ---------------- CONNECT TO BRAINAMPS ------------------%
        try
            % add any setup commands here
            if ~isfield(h,'BV'); 
                h.BV.delay = 0.001;
            end
            h.BV.status = 0;
        catch
            choice = questdlg('Brainvision not connected. Stop or continue?', ...
                    'Choice', ...
                    'Stop','Continue','Continue');
            switch choice
                case 'Stop'
                    error('Experiment stopped')
            end
        end
    case 'start'
        %% BV: start recording
        try
            %Priority(2);
            %NetStation('StartRecording');
            %pause(5)
            %NetStation('Synchronize');
            %NetStation('Event','BGIN');
            %pause(1);
        catch
            choice = questdlg('Brainvision not started. Stop or continue?', ...
                    'Choice', ...
                    'Stop','Continue','Continue');
            switch choice
                case 'Stop'
                    error('Experiment stopped')
            end
        end
    case 'stop'
        %% ---------------- STOP RECORDING ------------------%
        try
            %pause(1);
            %NetStation('Event','BEND');
            %pause(5);
            %Priority(0);
            %NetStation('StopRecording');
        catch
            msgbox('Brainvision is still recording','Message');
        end
    case 'mark'
        %NetStation('Event','STIM',stimtime,0.001,'TNUM',ns,'CNUM',h.Seq.condnum(h.i),'FNUM',h.Seq.signal(h.i),'BNUM',h.Seq.block(h.i));
        %--- EEG trigger on, via Labjack USB - for Brainamps ---%
        Error = ljud_AddRequest(h.ljHandle,LJ_ioPUT_DIGITAL_BIT,5,1,0,0);
        Error_Message(Error)
        
        % wait for a short delay
        Error = ljud_AddRequest(h.ljHandle,LJ_ioPUT_WAIT,5,h.BV.delay,0,0);
        Error_Message(Error)

        %--- EEG trigger off, via Labjack USB - for Brainamps ---%
        Error = ljud_AddRequest(h.ljHandle,LJ_ioPUT_DIGITAL_BIT,5,0,0,0);
        Error_Message(Error)
        
        % send command
        Error = ljud_GoOne(h.ljHandle);
        Error_Message(Error)
end