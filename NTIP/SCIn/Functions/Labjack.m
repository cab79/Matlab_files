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
                ME=0;
                while ME~=1 && t2<t1+60;
                    ME=connect_labjack(ME,t1,t2);
                    t2=GetSecs;
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
                
                % Configure LJTick-DAC, if needed
                try
                    Error = ljud_ePut(h.ljHandle, LJ_ioPUT_CONFIG, LJ_chTDAC_SCL_PIN_NUM,h.Settings.labjack_DACport,0); 
                    Error_Message(Error)
                    %Set DACA 
                    Error = ljud_ePut(h.ljHandle, LJ_ioTDAC_COMMUNICATION, LJ_chTDAC_UPDATE_DACA, 0, 0); 
                    Error_Message(Error)
                end

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

function ME=connect_labjack(ME,t1,t2)

try
    disp('loading labjack driver')
    ljud_LoadDriver; % Loads LabJack UD Function Library
    ME=1;
catch ME
    disp(['failed attempt'])
    pause(1)
end