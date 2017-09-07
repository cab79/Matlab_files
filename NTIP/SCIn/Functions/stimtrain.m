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
end