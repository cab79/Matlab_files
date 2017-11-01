function h = recordEEG(h,opt)

switch opt
    case 'connect'
        try
            if strcmp(h.Settings.record_EEG,'NS')
                %% ---------------- CONNECT TO NETSTATION ------------------%
                nshost = '10.0.0.42';
                nsport = 55513;

                if ~isfield(h,'NS')
                    h.NS=struc;
                end

                if ~isfield(h.NS,'status') && ...
                        exist('nshost','var') && ~isempty(nshost) && ...
                        exist('nsport','var') && nsport ~= 0
                    fprintf('Connecting to Net Station.\n');
                    [h.NS.status, nserror] = NetStation('Connect',nshost,nsport);
                    if h.NS.status ~= 0
                        error('Could not connect to NetStation host %s:%d.\n%s\n', ...
                            nshost, nsport, nserror);
                    end
                end
                NetStation('Synchronize');
            elseif strcmp(h.Settings.record_EEG,'BV')
                %% ---------------- CONNECT TO BRAINVISION ------------------%
                if ~isfield(h,'BV'); 
                    h.BV.delay = 0.001;
                end
                h.BV.status = 0;
                
            elseif strcmp(h.Settings.record_EEG,'serial')
                opt = 'EEG';
                open_serial(h,opt);
            end
            
            choice = questdlg('Using EEG: Disable GUI buttons? (more accurate EEG markers)', ...
                    'Choice', ...
                    'Yes','No','Yes');
            switch choice
                case 'Yes'
                    h.disableGUI=1;
            end
        catch MExc
            choice = questdlg('EEG not connected. Stop or continue?', ...
                    'Choice', ...
                    'Stop','Continue','Continue');
            switch choice
                case 'Stop'
                    error('Experiment stopped')
            end
        end
    case 'start'
        %% NS: start recording
        try
            if strcmp(h.Settings.record_EEG,'NS')
                Priority(2);
                NetStation('StartRecording');
                pause(5)
                NetStation('Synchronize');
                NetStation('Event','BGIN');
                pause(1);
            elseif strcmp(h.Settings.record_EEG,'BV')
                % do nothing
            end
        catch
            choice = questdlg('EEG not started. Stop or continue?', ...
                    'Choice', ...
                    'Stop','Continue','Continue');
            switch choice
                case 'Stop'
                    error('Experiment stopped')
            end
        end
    case 'pause'
        NetStation('Event','PAUS'); % mark EEG
        pause(1)
    case 'resume'
        NetStation('Event','CONT'); % mark EEG
        pause(1)
    case 'stop'
        %% ---------------- STOP RECORDING ------------------%
        try
            if strcmp(h.Settings.record_EEG,'NS')
                pause(1);
                NetStation('Event','STOP');
                pause(5);
                Priority(0);
                NetStation('StopRecording');
            elseif strcmp(h.Settings.record_EEG,'BV')
                % nothing
            end
        catch
            msgbox('EEG is still recording','Message');
        end
    case 'mark'
        if strcmp(h.Settings.record_EEG,'NS')
            NetStation('Event','STIM',stimtime,0.001,'TNUM',ns,'CNUM',h.Seq.condnum(h.i),'FNUM',h.Seq.signal(h.i),'BNUM',h.Seq.block(h.i));
            
        elseif strcmp(h.Settings.record_EEG,'BV')
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
        elseif strcmp(h.Settings.record_EEG,'serial')
            try
                global spt1
                TriggerNum = h.Seq.signal(h.i);
                fprintf(spt1,num2str(32+TriggerNum));
            catch
                disp('EEG not connected');
            end
        end
end