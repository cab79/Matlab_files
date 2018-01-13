function D=gplotprepare_eeglabdata_from_spm(S,D)
%% load EEGLAB data using filenames from an SPM file and reformat for plotting topography
% requires P containing fields:
    % P.spm_path
    % P.eeglab_path
    % P.eventtypes
    % P.savename
    % P.st_string
    % P.en_string

% load SPM
load(fullfile(D.spm_path,'SPM.mat'));
imglist = SPM.xY.P; % Subject-condition image list
fnames={};
for i = 1:length(imglist);
    st=strfind(imglist{i},S.st_string);
    en=strfind(imglist{i},S.en_string);
    fnames{i,1} =imglist{i}(st+length(S.st_string):en-1);
end
[fnames,~,fni] = unique(fnames,'stable');
D.DAT=cell(length(fnames),length(S.eventtypes));
for f = 1:length(fnames)
    EEGall=pop_loadset([fnames{f} '.set'],S.eeglab_path);
    for e = 1:length(S.eventtypes)
        if any(strcmp({EEGall.event.type},S.eventtypes{e}))
           EEG = pop_selectevent(EEGall,'type',S.eventtypes{e});
           D.DAT{f,e} = mean(EEG.data,3);
        end
    end
end
D.DAT=reshape(D.DAT,[],1);
D.EEGtimes = EEGall.times;
D.chanlocs=EEGall.chanlocs;
save(fullfile(S.eeglab_path,S.ERPsavename),'D');
