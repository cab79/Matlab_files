function handles = AuditoryAWE(handles,opt)

switch opt
    case 'start'
        play(handles.Seq)
        
    case 'pause'
        pause(handles.Seq);
    
    case 'resume'
        resume(handles.Seq);
    
    case 'stop'
        stop(handles.Seq);
end
