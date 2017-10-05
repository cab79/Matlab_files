function h = Labjack(h,opt)

switch opt
    case 'connect'
        %% ---------------- INITIALISING LabJack-----------------------------%
        try 
            if ~isfield(h,'ljHandle')

                fprintf('Connecting to LabJack.\n');

                ljud_LoadDriver; % Loads LabJack UD Function Library
                ljud_Constants; % Loads LabJack UD constant file
                [Error h.ljHandle] = ljud_OpenLabJack(LJ_dtU3,LJ_ctUSB,'1',1); % Returns ljHandle for open LabJack
                Error_Message(Error)

                %Start by using the pin_configuration_reset IOType so that all
                %pin assignments are in the factory default condition.
                Error = ljud_ePut (h.ljHandle, LJ_ioPIN_CONFIGURATION_RESET, 0, 0, 0);
                Error_Message(Error) % Checks for errors and displays them if they occur

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