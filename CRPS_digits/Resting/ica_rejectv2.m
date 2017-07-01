clear all
filepath = 'W:\Data\CRPS_resting\EEG';
%typerange = {'RHAN','LHAN','BELB','RELX'};
typerange = {'RHAN_subcomp','LHAN_subcomp','BELB_subcomp','RELX_subcomp'};
files = dir(fullfile(filepath,'*100Hz.Exp3.set'));
%[ALLEEG EEG CURRENTSET ALLCOM] = eeglab; % start EEGLAB under Matlab 

fnum = 16;

for f = fnum

    filename = files(f).name;
    [pth nme ext] = fileparts(filename);
    for e = 1:length(typerange)
        CURRENTSET = e;
        clear EEG
        EEG = pop_loadset('filename',[nme '_' typerange{e} ext],'filepath',filepath);
        EEG = pop_rejepoch(EEG, [], 0);
        ALLEEG(e) = EEG;
        if ~exist('rcomp','var'); 
            EEG = pop_selectcomps(EEG,1:35); 
            pause;
        end
        rcomp = EEG.reject.gcompreject;
        if ~exist('rejcomp','var'); rejcomp = rcomp; end
        rejcomp = (rcomp | rejcomp);
        clear EEG 
    end
    
    for e = 1:length(typerange)
        CURRENTSET = e;
        EEG = ALLEEG(e);
        EEG.reject.gcompreject = rejcomp;
        EEG = pop_subcomp(EEG,find(rejcomp),0);
        ALLEEG(e) = EEG;
        clear EEG
    end
    
    for e = 1:length(typerange)
        EEG = ALLEEG(e);
        CURRENTSET = e;
        assignin('base', 'EEG', EEG);
        evalin('base','detectArtefactsFinished=0;');
        
        cmd = [ ...
            '[tmprej tmprejE] = eegplot2trial( TMPREJ,' num2str(EEG.pnts) ',' num2str(EEG.trials) ');' ...
            'EEG = pop_rejepoch(EEG, tmprej, 1);' ...
            'detectArtefactsFinished=1;' ...
            ] ;

            %Draw the data.
            eegplot(EEG.data,...
                'srate',EEG.srate,...
                'events',EEG.event,...
                'command',cmd,...
                'butlabel','Reject');

            %Wait until the user has finished reviewing.
            reviewFinished=0;
            while ~reviewFinished
                reviewFinished=evalin('base','detectArtefactsFinished');
                pause(0.01);
            end
            evalin('base','clear detectArtefactsFinished');

            EEG=evalin('base','EEG');
        
        %pop_saveset(EEG,'filename',[nme '_' typerange{e} '_subcomp' ext],'filepath',filepath);
        pop_saveset(EEG,'filename',[nme '_' typerange{e} ext],'filepath',filepath);
        clear EEG;
    end
    
    ALLEEG = pop_delset(ALLEEG, [1:4]);
    clear EEG rejcomp rcomp ALLEEG
    
end

