function h = footpedal(h,opt)

h.buttonopt = {'1','2','3'};

switch opt
    case 'select'
        h.pedal = inputdlg('Foot pedal (1,2 or 3; leave blank if no foot pedal required): ');
    case 'start'
        disp('The experiment will start when after pressing the footswitch');
        keyIsDown=0;
        pause(1)
        while keyIsDown==0
            [keyIsDown,secs, keyCode] = KbCheck;
        end
end