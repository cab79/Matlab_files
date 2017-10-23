function h = blockstart(h,opt)
% command options to start each block

switch opt
    case 'msgbox'
        h.go = msgbox('Press OK to start the block.','Message');
        uiwait(h.go);
        
    case 'labjack'
        % RECEIVE SCANNER TRIGGER???
        choice = questdlg('Wait for scanner trigger(s) before starting?', ...
                'Choice', ...
                'Yes','No','Yes');
        switch choice
            case 'Yes'
                ljud_LoadDriver; % Loads LabJack UD Function Library
                ljud_Constants; % Loads LabJack UD constant file
                scantrig=[];
                lj_chan = 4;

                i=0;
                num_trig=0;
                h.trigtime=[];
                while num_trig<h.Settings.num_scanner_trig
                    [ljError] = ljud_ePut(h.ljHandle,LJ_ioPUT_DIGITAL_BIT,lj_chan,0,0);
                    scantrig=0;
                    %t1=GetSecs;
                    i=i+1;
                    disp(['Waiting for scanner trigger ' num2str(i)])
                    while scantrig==0
                        [ljError, scantrig] = ljud_eGet(h.ljHandle,LJ_ioGET_DIGITAL_BIT,lj_chan,0,0);

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
                    h.trigtime(i) = GetSecs;
                    num_trig = num_trig+1;
                    disp(['Triggered ' num2str(num_trig)]);
                    pause(h.Settings.waittime_scanner_trig);
                end
                disp('Triggers complete');
            case 'No'
                return
        end
        
    case 'buttonpress'
        opt = 'start';
        h = buttonpress(h,opt);
        
    
    case 'scannertrig'
        opt = 'scannertrig';
        h = buttonpress(h,opt);
        
    case 'audio'
        disp('PLAYING AUDIO')
        h = audio(h);

%pause(1);

end