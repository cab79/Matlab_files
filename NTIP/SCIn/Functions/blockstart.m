function h = blockstart(h,opt)
% command options to start each block

switch opt
    case 'msgbox'
        h.go = msgbox('Press OK to start the block.','Message');
        uiwait(h.go);
        
    case 'labjack'
        % RECEIVE SCANNER TRIGGER???
        disp('Waiting for scanner trigger');
        scantrig=[];
        lj_chan = 4;

        [ljError] = ljud_ePut (ljHandle,LJ_ioPUT_DIGITAL_PORT,lj_chan,0,0);
        while isempty(scantrig)
            [ljError, scantrig] = ljud_eGet (ljHandle,LJ_ioGET_DIGITAL_BIT_STATE,lj_chan,0,0);
        end
        Error_Message(ljError)
        disp('Triggered');
        
    case 'buttonpress'
        opt = 'start';
        h = buttonpress(h,opt);
        
    case 'audio'
        h = audio(h);

%pause(1);

end