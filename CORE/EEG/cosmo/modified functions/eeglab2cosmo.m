function cos = eeglab2cosmo(EEG,timewin,eventType)
% create cosmo data struct

if istruct(EEG)
    lats= dsearchn(EEG.times',timewin')';
    latidx=lats(1):lats(2);

    data = permute(EEG.data(:,latidx,:),[3 1 2]);
    cos.samples = reshape(data,size(data,1),[]);
    cos.fa.chan = repmat(1:EEG.nbchan,1,length(latidx));
    cos.fa.time = reshape(repmat(1:length(latidx),EEG.nbchan,1),1,[]);
    cos.a.fdim.labels = {'chan';'time'};
    cos.a.fdim.values = {{EEG.chanlocs.labels};EEG.times(latidx)/1000};
    cos.a.meeg.samples_field = 'trial';
    if isempty(eventType)
        conds = get_markers(EEG);
    else
        conds = eventType;
    end
    cos.sa.trialinfo(:,1) = conds';
    conduni = unique(conds);
    for c = 1:length(conduni)
        if isnumeric(conduni(c))
            condind = find(conds==conduni(c));
        else
            condind = strcmp(conds,conduni(c));
        end
        cos.sa.trialinfo(condind,2) = 1:length(condind)';
    end
    cos.labels=num2cell(cos.sa.targets);

elseif isnumeric(EEG)
    data=EEG;

    n_samples = size(data,2);
    n_chan = size(data,1);
    data = permute(data,[3 1 2]); % trials, chans/comps, samples
    cos.samples = reshape(data,size(data,1),[]); % trials, chans/comps*samples
    cos.fa.chan = repmat(1:n_chan,1,n_samples);
    cos.fa.time = reshape(repmat(1:n_samples,n_chan,1),1,[]);
    cos.a.fdim.labels = {'chan';'time'};
    cos.a.fdim.values = {cellstr(1:n_chan);timewin/1000};
    cos.a.meeg.samples_field = 'trial';
    cos.sa.trialinfo(:,1) = eventType';
    conduni = unique(cos.sa.trialinfo(:,1));
    for cn = 1:length(conduni)
        if isnumeric(conduni(cn))
            condind = find(cos.sa.trialinfo(:,1)==conduni(cn));
        else
            condind = strcmp(cos.sa.trialinfo(:,1),conduni(cn));
        end
        cos.sa.trialinfo(condind,2) = [1:length(condind)]';
    end
    cos.labels=num2cell(cos.sa.targets);
end