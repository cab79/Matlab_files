function h = D188(h,opt)

switch opt
    case 'setup'
        if ~isfield(h,'D188'); 
            h.d188.delay=1;  % delay between commands to D188 in ms
        end
    case 'outchan'
        % set output channel
        Error = ljud_AddRequest(h.ljHandle,LJ_ioPUT_DIGITAL_PORT,8,h.D188.chan,4,0);
        Error_Message(Error)
        
        % wait for a short delay
        Error = ljud_AddRequest(h.ljHandle,LJ_ioPUT_WAIT,4,h.D188.delay,0,0);
        Error_Message(Error)
        
        % send command
        Error = ljud_GoOne(h.ljHandle);
        Error_Message(Error)
end