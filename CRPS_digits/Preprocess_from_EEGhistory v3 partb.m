clear all
if isunix
    filepath = '/scratch/cb802/Data';
    run('/scratch/cb802/Matlab_files/CRPS_digits/loadsubj.m');
else
    filepath = 'W:\Data';
    run('W:\Matlab_files\CRPS_digits\loadsubj.m');
end
ana_path1 = fullfile(filepath,'CRPS_Digit_Perception_exp1');
ana_path2 = fullfile(filepath,'CRPS_Digit_Perception');
raw_path = fullfile(filepath,'CRPS_raw');
cd(raw_path);
files = dir('*ICA.set');

for f = 1:length(files)
    fname = files(f).name;
    [pth nme ext] = fileparts(fname);
    C = strsplit(nme,'_');
    EEG = pop_loadset('filename',fname,'filepath',raw_path);
    tname = [C{1} '_' C{2} '_' C{3} '_trials.mat'];
    load(tname);
    EEG = pop_subcomp( EEG, [], 0);
    EEG = pop_reref( EEG, []);
    
    if findstr(nme,'Exp1')
        if findstr(nme,'left')
            selectfnum =1:5;
            part1analysis;
        elseif findstr(nme,'right')
            selectfnum =6:10;
            part1analysis;
        end
    elseif findstr(nme,'Exp2')
        if findstr(nme,'left')
            selectfnum =1:5;
            part2analysis;
        elseif findstr(nme,'right')
            selectfnum =6:10;
            part2analysis;
        end
    end
   
    sname = [C{1} '_' C{2} '_' C{3} '.set'];
    EEG = pop_saveset(EEG,'filename',sname,'filepath',raw_path);
end