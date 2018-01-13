function h = Labjack(h,opt)

switch opt
    case 'connect'
        %% ---------------- INITIALISING LabJack-----------------------------%
        try 
            if ~isfield(h,'ljHandle')

                fprintf('Connecting to LabJack.\n');

                t1=GetSecs;
                t2=t1;
                att=0;
                while t2<t1+60;
                    try
                        att=att+1;
                        disp(['loading labjack driver, attempt ' num2str(att)])
                        t2=GetSecs;
                        ljud_LoadDriver; % Loads LabJack UD Function Library
                        break
                    catch
                        disp(['failed attempt'])
                        pause(0.5)
                    end
                end
                ljud_Constants; % Loads LabJack UD constant file
                [Error h.ljHandle] = ljud_OpenLabJack(LJ_dtU3,LJ_ctUSB,'1',1); % Returns ljHandle for open LabJack
                Error_Message(Error)

                %Start by using the pin_configuration_reset IOType so that all
                %pin assignments are in the factory default condition.
                Error = ljud_ePut (h.ljHandle, LJ_ioPIN_CONFIGURATION_RESET, 0, 0, 0);
                Error_Message(Error) % Checks for errors and displays them if they occur

                load LJfreqtable
                h.LJfreqtable = LJfreqtable;
            end
        catch
            choice = questdlg('Labjack not connected. Stop or continue?', ...
                    'Choice', ...
                    'Stop','Continue','Continue');
            switch choice
                case 'Stop'
                    error('Experiment stopped')
            end
        end
end