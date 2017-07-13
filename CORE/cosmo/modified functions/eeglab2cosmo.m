function cos = eeglab2cosmo(EEG,timewin)

lats= dsearchn(EEG.times',timewin')';
latidx=lats(1):lats(2);

data = permute(EEG.data(:,latidx,:),[3 1 2]);
cos.samples = reshape(data,size(data,1),[]);
cos.fa.chan = repmat(1:EEG.nbchan,1,length(latidx));
cos.fa.time = reshape(repmat(1:length(latidx),EEG.nbchan,1),1,[]);
cos.a.fdim.labels = {'chan';'time'};
cos.a.fdim.values = {{EEG.chanlocs.labels};EEG.times(latidx)/1000};
cos.a.meeg.samples_field = 'trial';
conds = get_markers(EEG);
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