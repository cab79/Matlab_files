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
        %keyIsDown=0;
        %while keyIsDown==0
        %    keyIsDown = KbCheck;
        %end
        try % if there is already a queue
            KbQueueWait; 
            KbQueueFlush;
        catch
            KbQueueCreate;	
            KbQueueStart;
            KbQueueWait; 
            KbQueueRelease;
        end
        disp('Button pressed: Starting experiment');
        set(GUIh.info, 'String', 'Running sequence...');
        
    case 'scannertrig'
        scantrig=[];
        GUIh = guidata(h.GUIhname);

        i=0;
        num_trig=0;
        h.trigtime=[];
        while num_trig<h.Settings.num_scanner_trig
            scantrig=0;
            %t1=GetSecs;
            i=i+1;
            disp(['Waiting for scanner trigger ' num2str(i)])
            while scantrig==0 
                keyIsDown=0;
                [keyIsDown, presstime, keyCode, deltaSecs] = KbCheck;
                if keyIsDown
                    btnstr = KbName(keyCode);
                    if ismember(btnstr,h.Settings.triggeropt)
                        scantrig =1;
                    end
                end
                pause(0.001)

                %% for testing only:
                %pause(1)
                %t2=GetSecs;
                %if t2-t1>=5
                %    [ljError] = ljud_ePut(h.ljHandle,LJ_ioPUT_DIGITAL_BIT,5,1,0);
                %end
                %% 
            end
            %% for testing only:
            %[ljError] = ljud_ePut(h.ljHandle,LJ_ioPUT_DIGITAL_BIT,5,0,0);
            %Error_Message(ljError)
            %%
            %h.trigtime(i) = presstime;
            
            h.out.presstrial = [h.out.presstrial 0];
            h.out.pressbutton = [h.out.pressbutton {btnstr}];
            h.out.presstime = [h.out.presstime presstime];
            h.out.presstimedelta = [h.out.presstimedelta deltaSecs];
            h.out.RT = [h.out.RT 0];
            
            num_trig = num_trig+1;
            disp(['Triggered ' num2str(num_trig)]);
            pause(h.Settings.waittime_scanner_trig);
        end
        disp('Triggers complete');
end