function h = Netstation(h,opt)

switch opt
    case 'connect'
        %% ---------------- CONNECT TO NETSTATION ------------------%
        try
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
        catch
            choice = questdlg('Netstation not connected. Stop or continue?', ...
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
            Priority(2);
            NetStation('StartRecording');
            pause(5)
            NetStation('Synchronize');
            NetStation('Event','BGIN');
            pause(1);
        catch
            choice = questdlg('Netstation not started. Stop or continue?', ...
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
            pause(1);
            NetStation('Event','STOP');
            pause(5);
            Priority(0);
            NetStation('StopRecording');
        catch
            msgbox('Netstation is still recording','Message');
        end
    case 'mark'
        NetStation('Event','STIM',stimtime,0.001,'TNUM',ns,'CNUM',h.Seq.condnum(h.i),'FNUM',h.Seq.signal(h.i),'BNUM',h.Seq.block(h.i));
end