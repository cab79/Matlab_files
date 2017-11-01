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

if option==2
    %Extract only S1 (standard)
    EEG = pop_loadset('filename', filename,'filepath', filepath);
    EEG = eeg_checkset( EEG );    
    EEG = pop_selectevent( EEG, 'type',{'S  1'},'deleteevents','on','deleteepochs','on','invertepochs','off');
    %Save new set
    sname = [C{1} '_' C{2} '_S1.set'];
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset(EEG,'filename',sname,'filepath',anapath);
    
    %Extract only S2
    EEG = pop_loadset('filename', filename,'filepath', filepath);
    EEG = eeg_checkset( EEG );    
    EEG = pop_selectevent( EEG, 'type',{'S  2'},'deleteevents','on','deleteepochs','on','invertepochs','off');
    %Save new set
    sname = [C{1} '_' C{2} '_S2.set'];
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset(EEG,'filename',sname,'filepath',anapath);
    
    %Extract only S3
    EEG = pop_loadset('filename', filename,'filepath', filepath);
    EEG = eeg_checkset( EEG );    
    EEG = pop_selectevent( EEG, 'type',{'S  3'},'deleteevents','on','deleteepochs','on','invertepochs','off');
    %Save new set
    sname = [C{1} '_' C{2} '_S3.set'];
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset(EEG,'filename',sname,'filepath',anapath);    
    %Plot channel ERP in scalp array
    
    %Extract only S4
    EEG = pop_loadset('filename', filename,'filepath', filepath);
    EEG = eeg_checkset( EEG );    
    EEG = pop_selectevent( EEG, 'type',{'S  4'},'deleteevents','on','deleteepochs','on','invertepochs','off');
    %Save new set
    sname = [C{1} '_' C{2} '_S4.set'];
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset(EEG,'filename',sname,'filepath',anapath); 
    
    %Extract only S5
    EEG = pop_loadset('filename', filename,'filepath', filepath);
    EEG = eeg_checkset( EEG );    
    EEG = pop_selectevent( EEG, 'type',{'S  5'},'deleteevents','on','deleteepochs','on','invertepochs','off');
    %Save new set
    sname = [C{1} '_' C{2} '_S5.set'];
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset(EEG,'filename',sname,'filepath',anapath); 
    
    %Extract only S6
    EEG = pop_loadset('filename', filename,'filepath', filepath);
    EEG = eeg_checkset( EEG );    
    EEG = pop_selectevent( EEG, 'type',{'S  6'},'deleteevents','on','deleteepochs','on','invertepochs','off');
    %Save new set
    sname = [C{1} '_' C{2} '_S6.set'];
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset(EEG,'filename',sname,'filepath',anapath);     
    
    %Extract only S7
    EEG = pop_loadset('filename', filename,'filepath', filepath);
    EEG = eeg_checkset( EEG );    
    EEG = pop_selectevent( EEG, 'type',{'S  7'},'deleteevents','on','deleteepochs','on','invertepochs','off');
    %Save new set
    sname = [C{1} '_' C{2} '_S7.set'];
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset(EEG,'filename',sname,'filepath',anapath);     
    
    %Extract only S8
    EEG = pop_loadset('filename', filename,'filepath', filepath);
    EEG = eeg_checkset( EEG );    
    EEG = pop_selectevent( EEG, 'type',{'S  8'},'deleteevents','on','deleteepochs','on','invertepochs','off');
    %Save new set
    sname = [C{1} '_' C{2} '_S8.set'];
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset(EEG,'filename',sname,'filepath',anapath);    
    
    %Extract only S9
    EEG = pop_loadset('filename', filename,'filepath', filepath);
    EEG = eeg_checkset( EEG );    
    EEG = pop_selectevent( EEG, 'type',{'S  9'},'deleteevents','on','deleteepochs','on','invertepochs','off');
    %Save new set
    sname = [C{1} '_' C{2} '_S9.set'];
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset(EEG,'filename',sname,'filepath',anapath);    
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
