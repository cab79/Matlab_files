% ERP EXTRACTION
% Compare the ERP elicited by 'standard' stimulus to that elicited by the 'deviant'
% stimulus. Two options depending on what you want to plot to compare:
%               - Option 1: compare 'standard' ERP to the combined average of the 'deviants' ERP's. 
%               - Option 2: compare 'standard' ERP to each individual 'deviant' ERP.

%File locations
filepath = 'Y:\Marie Shorrock\NTIP\Pilot_Tim_Auditory\Preprocessed';
filename = 'NTIP_TimAudioEnableGUI_cleaned.set';
anapath = 'Y:\Marie Shorrock\NTIP\Pilot_Tim_Auditory\Analysis\ERP Separation'; % folder to save analysed data

cd(filepath);
cd(anapath);

% Choose the option (1 or 2)
option = 1;

if option==1
    %Extract only S1 (standard)
    EEG = pop_loadset('filename', filename,'filepath', filepath);
    EEG = eeg_checkset( EEG );    
    EEG = pop_selectevent( EEG, 'type',{'S  1'},'deleteevents','on','deleteepochs','on','invertepochs','off');
    %Save new set
    sname = [C{1} '_' C{2} '_S1.set'];
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset(EEG,'filename',sname,'filepath',anapath);
    
    %Extract all deviants and remove S1
    EEG = pop_loadset('filename', filename,'filepath', filepath);
    EEG = eeg_checkset( EEG );
    EEG = pop_selectevent( EEG, 'type',{'S  2' 'S  3' 'S  4' 'S  5' 'S  6' 'S  7' 'S  8' 'S  9'},'deleteevents','on','deleteepochs','on','invertepochs','off');
    %Save new set
    sname = [C{1} '_' C{2} '_SAD.set'];
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset(EEG,'filename',sname,'filepath',anapath);
end


% EXAMPLE OF HOW TO MAKE THIS INTO A LOOP:
if option==2
    stim_num = {'1','2','3','4','5','6','7','8','9'};
    %load
    EEGallstim = pop_loadset('filename', filename,'filepath', filepath);
    EEGallstim = eeg_checkset( EEGallstim );  
    for sn = 1:length(stim_num)
        %Extract  
        EEG = pop_selectevent( EEGallstim, 'type',{['S  ' stim_num{sn}]},'deleteevents','on','deleteepochs','on','invertepochs','off');
        %Save new set
        sname = [C{1} '_' C{2} '_S' stim_num{sn} '.set'];
        EEG = eeg_checkset( EEG );
        EEG = pop_saveset(EEG,'filename',sname,'filepath',anapath);
    end
   
end

% EXAMPLE OF HOW TO PLOT WITHOUT SAVING SEPARATE FIRST FILES (saves harddrive space!):
if option==2
    stim_num = {'1','2','3','4','5','6','7','8','9'};
    
    %load
    EEGallstim = pop_loadset('filename', filename,'filepath', filepath);
    EEGallstim = eeg_checkset( EEGallstim );  
    
    for sn = 1:length(stim_num)
        %Extract  
        EEG = pop_selectevent( EEGallstim, 'type',{['S  ' stim_num{sn}]},'deleteevents','on','deleteepochs','on','invertepochs','off');
        %Plot
        figure; pop_plottopo(EEG, [1:62] , 'files', 0, 'ydir', 1, 'limits', [-200 299 -2.25 2.25]);
        figure; pop_topoplot(EEG,1, 299,'',[1 1] ,0,'electrodes','on', 'maplimits', [-0.7, 0.7]);
        eeglab redraw; 
    end
   
end

%% FOR PLOTTING
%%Find all the files that have been saved in a different file
filepath = 'Y:\Marie Shorrock\NTIP\Pilot_Tim_Auditory\Analysis\ERP Separation';
cd(filepath);
files = dir('NTIP_TimAudioEnableGUI_S*.set')
files_ana = 1:length(files);

%Plot channel ERP's in scalp array
for f = files_ana
    [pth nme ext] = fileparts(files(f).name); % extract parts of filename for use later, e.g. saving
    C = strsplit(nme,'_');
    EEG = pop_loadset('filename',files(f).name,'filepath',filepath);    
    EEG = eeg_checkset( EEG );    
    figure; pop_plottopo(EEG, [1:62] , 'files', 0, 'ydir', 1, 'limits', [-200 299 -2.25 2.25]);
    eeglab redraw; 
end     
