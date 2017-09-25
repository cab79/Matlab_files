function h = buttonpress(h,opt)


switch opt
    case 'setup'
        %% initialise GUI button info
        %GUIh = guidata(h.GUIhname);
        %set(GUIh.buttonpressed, 'String', 'Waiting');
        %set(GUIh.buttontime, 'String', 0);
        %guidata(h.GUIhname,GUIh)
    case 'start'
        GUIh = guidata(h.GUIhname);
        disp('The experiment will start after pressing the button');
        set(GUIh.info, 'String', 'Press button to start');
        pause(0.001)
        %% get GUI data
        %pressed=0;
        %while pressed==0
            %GUIh = guidata(h.GUIhname);
        %    pause(0.01)
        %    if ~strcmp(get(GUIh.buttonpressed,'String'),'Waiting')
        %        disp('Button pressed: Starting experiment');
        %        pressed=1;
        %    end
        %end
        
        %% to do the same thing using Psychotoolbox rather than GUI:
        keyIsDown=0;
        while keyIsDown==0
            keyIsDown = KbCheck;
        end
        disp('Button pressed: Starting experiment');
        set(GUIh.info, 'String', 'Running sequence...');
end