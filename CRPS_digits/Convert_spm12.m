clear all

filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception_exp1\alltrials\';
run('C:\Data\Matlab\Matlab_files\CRPS_digits\loadsubj.m');
cd(filepath)
files = dir('*t.set');
flip=1;

for f = 1:length(files)
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
    EEG.dataset = [nme '_change.set'];
    %prefix = 'spm8_dc_';
    prefix = 'spm12_';
    if flip==1 && (~isempty(strfind(nme, 'right')) || ~isempty(strfind(nme, 'Right'))); 
        EEG = flipchan(EEG); 
        prefix = 'spm12_flip_';
        %prefix = 'spm8_flip_';
    end
    
    EEG.outfile = [prefix spm_str_manip(EEG.dataset,'tr')];
    
    fevents = zeros(1,size(EEG.data,3));
        %j=0;
        for i = 1:size(EEG.data,3)
            if sum(strcmp(EEG.epoch(i).eventtype(1,:),'STIM'))>0
                ind = find(strcmp(EEG.epoch(i).eventtype(1,:),'STIM'));
                findfnum=[];
                %findcnum=[];
                if iscell(EEG.epoch(i).eventinit_index)
                    findfnum = find(strcmp('FNUM',EEG.epoch(i).eventcodes{1,ind(1)}(:,1)));
                    %findcnum = find(strcmp('CNUM',EEG.epoch(i).eventcodes{1,ind(1)}(:,1)));
                    fevents(1,i) = EEG.epoch(i).eventcodes{1,ind(1)}{findfnum(end),2};
                %    cevents(1,i) = EEG.epoch(i).eventcodes{1,ind(1)}{findcnum(end),2};
                else 
                    findfnum = find(strcmp('FNUM',EEG.epoch(i).eventcodes(:,ind(1))));
                    %findcnum = find(strcmp('CNUM',EEG.epoch(i).eventcodes(:,ind(1)))); 
                    fevents(1,i) = EEG.epoch(i).eventcodes{findfnum(end),2};
                    %cevents(1,i) = EEG.epoch(i).eventcodes{findcnum(end),2};
                end
            end
        end
    
    EEG.conditionlabels = cellfun(@num2str, num2cell(fevents), 'UniformOutput', false);
    EEG.mode ='epoched';
    EEG.timewin = [-200 800];
    spm_eeg_convert(EEG);
    %end
    clear EEG;
end
matlabmail
%x=1;save('x.mat','x');