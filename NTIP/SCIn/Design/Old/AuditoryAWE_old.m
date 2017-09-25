function h = AuditoryAWE(h,opt)

switch opt
    case 'start'
        play(h.Seq.signal)
        
    case 'pause'
        pause(h.Seq.signal);
    
    case 'resume'
        resume(h.Seq.signal);
    
    case 'stop'
        stop(h.Seq.signal);
end
