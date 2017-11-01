%LOAD DATA FOR CHANLOCS

%Original dataset location
origpath = 'Y:\Marie Shorrock\NTIP\Pilot_Tim_Auditory\Set'; % unprocessed data

% FIND THE DATA FILES
cd(origpath);
files = dir('NTIP_TimEyesClosed_orig.set');

files_ana = 1:length(files);
trials_ana = 1; fname_ext = '';

for f = files_ana
    
    orig_file = files(f).name;
    [pth nme ext] = fileparts(orig_file); 
    C = strsplit(nme,'_');
    
    % LOAD DATA
    EEG = pop_loadset('filename',orig_file,'filepath',origpath);
    
    chanlocs = EEG.chanlocs;
    save chanlocs chanlocs
end
