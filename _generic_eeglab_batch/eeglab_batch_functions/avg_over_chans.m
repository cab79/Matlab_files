function S=avg_over_chans(S)

S.func = 'ca';

% select channel groups
if ischar(S.(S.func).chanaverage.groups{1})
    % use chanlocs file
    chanlocs = load(S.(S.func).chanaverage.groups{1});
    field = getfields(chanlocs)
    chanlocs = chanlocs.(field);
    S=select_chans(S);
    chanlocs = chanlocs(S.(S.func).inclchan);
    uchangroups = unique(chanlocs);
    for g = 1:length(uchangroups)
        changroups{1,g} = find(chanlocs==uchangroups(g));
    end
elseif isnumeric(S.(S.func).chanaverage.groups{1})
    % or use channel indices from script
    changroups = S.(S.func).chanaverage.groups;
    S.(S.func).inclchan = [changroups{:}];
end

% get path to data
switch S.(S.func).select.datatype
    case 'ERP'
        S.path.file = S.path.erp;
    case 'TF'
        S.path.file = S.path.tf;
    case 'Freq'
        S.path.file = S.path.freq;
end
   

% GET FILE LIST
S = getfilelist(S);

% load data
for f = 1:length(S.(S.func).filelist)
    file = S.(S.func).filelist{f};
    if strfind(S.(S.func).fname.ext{:},'mat')
        % either ERP, TF or Freq data
        dat = load(fullfile(S.path.file,file));
        switch S.(S.func).select.datatype
            case 'ERP'
                if isfield(dat,'gadata')
                    dat.tldata = {dat.gadata};
                end
                for ev = 1:length(dat.tldata)
                    data = dat.tldata{ev}.avg;
                    for g = 1:length(changroups)
                        av_data{g} = squeeze(mean(data(changroups{g},:,:),1));
                    end
                    % create new data struct
                    dat.tldata{ev}.powspctrm=cat(ndims(data),av_data{:});
                    dat.tldata{ev}.powspctrm=permute(dat.tldata{ev}.powspctrm,[1 ndims(data) 2]);
                    label = {dat.tldata{ev}.label{S.(S.func).inclchan}};
                    dat.tldata{ev}.label = {};
                    for g = 1:length(changroups)
                        dat.tldata{ev}.label{g} = label(changroups{g});
                    end
                    tldata{ev} = dat.tldata{ev};
                end
                % save
                [~,nme,~] = fileparts(file);
                sname = [nme '_chanavg.mat'];
                if isfield(dat,'gadata')
                    gadata = tldata; % grandavg data
                    save(fullfile(S.path.file,sname),'gadata');
                else
                    save(fullfile(S.path.file,sname),'tldata');
                end
            case {'TF','Freq'}
                if isfield(dat,'gadata')
                    dat.fdata = {dat.gadata};
                end
                for ev = 1:length(dat.fdata)
                    data = dat.fdata{ev}.powspctrm;
                    if strfind(dat.fdata{ev}.dimord,'rpt')
                        for g = 1:length(changroups)
                            av_data{g} = squeeze(mean(data(:,changroups{g},:),2));
                        end
                        % create new data struct
                        dat.fdata{ev}.powspctrm=cat(ndims(data),av_data{:});
                        dat.fdata{ev}.powspctrm=permute(dat.fdata{ev}.powspctrm,[1 ndims(data) 2]);
                    else;
                        for g = 1:length(changroups)
                            av_data{g} = squeeze(mean(data(changroups{g},:,:),1));
                        end
                        % create new data struct
                        dat.fdata{ev}.powspctrm=cat(ndims(data),av_data{:});
                        dat.fdata{ev}.powspctrm=permute(dat.fdata{ev}.powspctrm,[ndims(data)+1 1 2]);
                    end;
                    
                    label = {dat.fdata{ev}.label{S.(S.func).inclchan}};
                    dat.fdata{ev}.label = {};
                    for g = 1:length(changroups)
                        dat.fdata{ev}.label{g} = label(changroups{g});
                    end
                    fdata{ev} = dat.fdata{ev};
                    if isfield(fdata{ev},'madata')
                        madata=[];
                        for g = 1:length(changroups)
                            madata(g,:,:) = squeeze(mean(fdata{ev}.madata(changroups{g},:,:),1));
                        end
                        fdata{ev}.madata = madata;
                    end
                end
                % save
                [~,nme,~] = fileparts(file);
                sname = [nme '_chanavg.mat'];
                if isfield(dat,'gadata')
                    gadata = fdata; % grandavg data
                    save(fullfile(S.path.file,sname),'gadata');
                else
                    save(fullfile(S.path.file,sname),'fdata');
                end
                
        end
    elseif strfind(S.(S.func).fname.ext,'set')
        % single-trial EEGLAB
        EEG = pop_loadset('filename',file,'filepath',S.path.file);
        data = EEG.data;
        for g = 1:length(changroups)
            av_data{g} = squeeze(mean(data(changroups{g},:,:),1));
        end
        % create new data struct
        EEG.data=av_data;
        label = EEG.chanlocs(S.(S.func).inclchan);
        EEG.chanlocs = {};
        for g = 1:length(changroups)
            EEG.chanlocs{g} = label(changroups{g});
        end
        % save
        [~,nme,~] = fileparts(file);
        sname = [nme '_chanavg.set'];
        pop_saveset(EEG,'filename',sname,'filepath',S.path.file);
    end
end

