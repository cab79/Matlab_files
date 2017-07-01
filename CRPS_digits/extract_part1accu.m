clear all

if isunix
    filepath = '/scratch/cb802/Data/CRPS_raw/Raw';
    run('/scratch/cb802/Matlab_files/CRPS_digits/loadsubj.m');
else
    filepath = 'W:\Data\CRPS_raw\Raw';
    run('W:\Matlab_files\CRPS_digits\loadsubj.m');
end

cd(filepath)
files = dir('*orig.set');

for f = 1:length(files)
    EEG = pop_loadset(files(f).name,filepath);
    [pth nme ext] = fileparts(EEG.filename);
    [pth orig_datfile ext] = fileparts(EEG.filename)
    orig_datfile = [orig_datfile '.fdt'];
    EEG.datfile = [nme '.fdt'];
    if isfield(EEG,'dataset'); EEG.dataset = [nme ext];end
    if exist('fullfile(filepath, orig_datfile)','file');
        movefile(fullfile(filepath, orig_datfile),fullfile(filepath, EEG.datfile));
    end
    if ~isfield(EEG,'dataset')
        [pth nme ext] = fileparts(EEG.datfile);
        EEG.dataset = [nme '.set'];
    end

    %calculates accuracy for each trial
    accudata = [];
    for e = 1:length(EEG.event)
        if strcmp(EEG.event(e).type,'STIM')
            thisfnum = EEG.event(e).codes{strcmp('FNUM',EEG.event(e).codes(:,1)),2};
            if e==length(EEG.event)
                thisrnum = thisfnum;
            else
                try
                    thisrnum = EEG.event(e+1).codes{strcmp('RNUM',EEG.event(e+1).codes(:,1)),2};
                catch
                    thisrnum = nan;
                end
            end
            accu = double(thisfnum == thisrnum);
            %if sum(strcmp('ACCU',EEG.event(e).codes(:,1))) == 0
            %    EEG.event(e).codes = cat(1,EEG.event(e).codes,{'ACCU' accu});
            %else
            %    EEG.event(e).codes{strcmp('ACCU',EEG.event(e).codes(:,1)),2} = accu;
            %end
            accudata = cat(1,accudata,[thisfnum thisrnum accu]);
        end
    end
    xlswrite(sprintf('%s_accu.xls',nme),accudata);
    clear EEG
end