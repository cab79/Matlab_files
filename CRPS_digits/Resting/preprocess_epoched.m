clear all
loadpaths
%filepath = 'W:\Data\CRPS_resting\EEG';
removechan
typerange = {'RHAN_subcomp','LHAN_subcomp','BELB_subcomp','RELX_subcomp'};
%typerangeS = {'RHAN','LHAN','BELB','RELX'};
%epoch_dur = 4;
files = dir(fullfile(filepath,'*100Hz.Exp3.set'));
%hicutoff = 45;

for f = 16%:length(files)
    
    filename = files(f).name;
    [pth nme ext] = fileparts(filename);
    
    for e = 1:length(typerange)
        efilename = [nme '_' typerange{e} ext];
        EEG = pop_loadset('filename',efilename,'filepath',filepath);
        ALLEEG(e) = EEG;
    end
    
    ALLEEG = pop_runica(ALLEEG,'icatype','jader','concatenate','on','options',{35});
   %ALLEEG = pop_runica(ALLEEG,'icatype','runica','concatenate','on');
    
    for e = 1:length(typerange)
        efilename = [nme '_' typerange{e} ext];
        pop_saveset(ALLEEG(e), efilename, filepath);
    end
    
    clear TOTEEG EEG ALLEEG
end

