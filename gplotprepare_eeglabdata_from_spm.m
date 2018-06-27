function E=gplotprepare_eeglabdata_from_spm(S,D)
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
fnamesdesign={};
for i = 1:length(imglist);
    st=strfind(imglist{i},S.st_string);
    en=strfind(imglist{i},S.en_string);
    fnamesdesign{i,1} =imglist{i}(st+length(S.st_string):en-1);
    
    %get imglist condition order
    [~,nme,~] = fileparts(imglist{i});
    C=strsplit(nme,'_');
    condlabelsdesign(i,1) = str2num(C{2});
end
condlabels = unique(condlabelsdesign,'stable');
condlabels = condlabels(ismember(condlabels,S.eventtypes));
[fnames,~,fni] = unique(fnamesdesign,'stable');
D.DAT=cell(length(imglist),1);
D.filenames=cell(length(imglist),2)
for f = 1:length(fnames)
    EEGall=pop_loadset([fnames{f} '.set'],S.eeglab_path);
    
    % obtain trial indices from EGI data and flip chans if needed
    if any(find(strcmp('STIM',EEGall.epoch(1).eventtype)))
        [conds, tnums, fnums, bnums] = get_markers(EEGall);
        if S.use_flipped
            EEGflip = flipchanEGI(EEGall);
            for i = S.flipcond
                trialind = find(conds==i);
                EEGall.data(:,:,trialind)= EEGflip.data(:,:,trialind);
            end
            clear EEGflip
        end
    end
    
    filecondind = find(fni==f);
    for e = 1:length(condlabels)
        
        if iscell(condlabels) && any(strcmp({EEGall.event.type},condlabels{e}))
           EEG = pop_selectevent(EEGall,'type',condlabels{e});
           D.DAT{filecondind(e),1} = mean(EEG.data,3);
        elseif exist('conds','var') % EGI
           EEG = pop_select(EEGall,'trial',find(conds==condlabels(e)));
           D.DAT{filecondind(e),1} = mean(EEG.data,3);
        end
        D.filenames(filecondind(e),:) = {fnames{f},condlabels(e)};
    end
end
D.DAT=reshape(D.DAT,[],1);
D.EEGtimes = EEGall.times;
D.chanlocs=EEGall.chanlocs;
E=D;
save(fullfile(S.eeglab_path,S.ERPsavename),'E');
