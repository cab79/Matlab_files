clear all
%if isunix
%    filepath = '/scratch/cb802/Data/CRPS_Digit_Perception';
%    run('/scratch/cb802/Matlab_files/CRPS_digits/loadsubj.m');
%else
%    filepath = 'W:\Data\CRPS_Digit_Perception';
%    run('W:\Matlab_files\CRPS_digits\loadsubj.m');
%end
if isunix
    filepath = '/scratch/cb802/Data/CRPS_Digit_Perception_exp1/alltrials';
    run('/scratch/cb802/Matlab_files/CRPS_digits/loadsubj.m');
else
    filepath = 'W:\Data\CRPS_Digit_Perception_exp1\alltrials';
    run('W:\Matlab_files\CRPS_digits\loadsubj.m');
end
cd(filepath)
files = dir('*t.set');
flip=1;

for f = 30:length(files)
    if any(strfind(files(f).name,'flip')) || any(strfind(files(f).name,'ICA')) || any(strfind(files(f).name,'orig'))
        continue;
    end;
    EEG = pop_loadset(files(f).name,filepath);
    EEG.filename = strrep(EEG.filename, '.left', '_left');
    EEG.filename = strrep(EEG.filename, '.Left', '_left');
    EEG.filename = strrep(EEG.filename, '.right', '_right');
    EEG.filename = strrep(EEG.filename, '.Right', '_right');
    EEG.filename = strrep(EEG.filename, '.flip', '_flip');
    EEG.filename = strrep(EEG.filename, '.Flip', '_flip');
    EEG.filename = strrep(EEG.filename, '.aff', '_aff');
    EEG.filename = strrep(EEG.filename, '.Aff', '_aff');
    EEG.filename = strrep(EEG.filename, '.Unaff', '_unaff');
    EEG.filename = strrep(EEG.filename, '.unaff', '_unaff');
    EEG.filename = strrep(EEG.filename, '_Left', '_left');
    EEG.filename = strrep(EEG.filename, '_Right', '_right');
    EEG.filename = strrep(EEG.filename, '_Flip', '_flip');
    EEG.filename = strrep(EEG.filename, '_Aff', '_aff');
    EEG.filename = strrep(EEG.filename, '_Unaff', '_unaff');
    EEG.filename = strrep(EEG.filename, '.Exp1', '_Exp1');
    [pth nme ext] = fileparts(EEG.filename);
    [pth orig_datfile ext] = fileparts(EEG.filename)
    orig_datfile = [orig_datfile '.fdt'];
    EEG.datfile = [nme '.fdt'];
    if isfield(EEG,'dataset'); EEG.dataset = [nme ext];end
    if exist('fullfile(filepath, orig_datfile)','file')
        movefile(fullfile(filepath, orig_datfile),fullfile(filepath, EEG.datfile));
    end
    [pth nme ext] = fileparts(EEG.datfile);
    EEG.dataset = [nme '.set'];
    %prefix = 'spm8_dc_';
    prefix = 'spm8_dc_';
    if flip==1 && (~isempty(strfind(nme, 'right')) || ~isempty(strfind(nme, 'Right'))); 
        EEG = flipchan(EEG); 
        prefix = 'spm8_dc_flip_';
        %prefix = 'spm8_flip_';
    end
    
    EEG.outfile = [prefix spm_str_manip(EEG.dataset,'tr')];
    spm_eeg_convert_ui(EEG,'justread');
    %end
    clear EEG;
end
matlabmail
%x=1;save('x.mat','x');