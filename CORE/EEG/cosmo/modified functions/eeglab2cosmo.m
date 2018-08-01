function cos = eeglab2cosmo(EEG,timewin,condlist)
% create cosmo data struct

if isstruct(EEG)
    lats= dsearchn(EEG.times',timewin')';
    latidx=lats(1):lats(2);

    data = permute(EEG.data(:,latidx,:),[3 1 2]);
    cos.samples = reshape(data,size(data,1),[]);
    cos.fa.chan = repmat(1:EEG.nbchan,1,length(latidx));
    cos.fa.time = reshape(repmat(1:length(latidx),EEG.nbchan,1),1,[]);
    cos.a.fdim.labels = {'chan';'time'};
    cos.a.fdim.values = {{EEG.chanlocs.labels};EEG.times(latidx)/1000};
    cos.a.meeg.samples_field = 'trial';
    if isempty(condlist)
        conds = get_markers(EEG);
    else
        conds = condlist;
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

elseif isnumeric(EEG)
    data=EEG;
    n_samples = size(data,2);
    n_chan = size(data,1);
    chan = 1:n_chan;
    chan_labels = cellfun(@num2str,num2cell(chan),'UniformOutput',false);
    chan_labels = strcat('E',chan_labels);

    data = permute(data,[3 1 2]); % trials, chans/comps, samples
    cos.samples = reshape(data,size(data,1),[]); % trials, chans/comps*samples
    cos.fa.chan = repmat(chan,1,n_samples);
    cos.fa.time = reshape(repmat(1:n_samples,n_chan,1),1,[]);
    cos.a.fdim.labels = {'chan';'time'};
    cos.a.fdim.values = {chan_labels;timewin/1000};
    cos.a.meeg.samples_field = 'trial';
    cos.sa.trialinfo(:,1) = condlist';
    conduni = unique(cos.sa.trialinfo(:,1));
    for cn = 1:length(conduni)
        if isnumeric(conduni(cn))
            condind = find(cos.sa.trialinfo(:,1)==conduni(cn));
        else
            condind = strcmp(cos.sa.trialinfo(:,1),conduni(cn));
        end
        cos.sa.trialinfo(condind,2) = [1:length(condind)]';
    end
end

% remove NaNs
rm = find(isnan(cos.sa.trialinfo(:,1)));
cos.sa.trialinfo(rm,:) = [];
cos.samples(rm,:) = [];