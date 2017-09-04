function h = buttonpress(h,opt)


switch opt
    case 'setup'
        h.button.select = inputdlg('Button to use for sujbect response; leave blank if no button press required): ');
        h.button.press2start = 1;
        h.button.opt = {'1','2','3'}; % range of allowable responses
    case 'start'
        if h.button.press2start
            disp('The experiment will start when after pressing the button');
            pressed=0;
            % get GUI data
            while pressed==0
                GUIh = guidata(h.GUIhname);
                if isfield(GUIh.button,'press')
                    GUIh = rmfield(GUIh.button,'press');
                    guidata(h.GUIhname,GUIh)
                    disp('Button pressed: Starting experiment');
                    pressed=1;
                end
            end
        end
        
        %% to do the same thing using Psychotoolbox rather than GUI:
        %keyIsDown=0;
        %pause(1)
        %while keyIsDown==0
        %    [keyIsDown,secs, keyCode] = KbCheck;
        %end
end