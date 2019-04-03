function S = flat_channel_reject(S,EEG)
% function maximises the amount of data remaining after removing flat
% channels (on certain trials)

% Planned future updates:
% add topographic plots showing locations of channels to be interpolated.

% OLD
% std threshold - less variance than this per trial will be rejected
% varthresh = S.prep.clean.flatchan.varthresh;

results = [];
idx={};

% inverse variance over data points on each trial and channel
invvar = 1./squeeze(std(EEG.data,[],2));

% order variance over trials: trialind is an index of which trials have the
% least variance
[~,trialind]=sort(sum(invvar,1),'descend');

% cycle through every combination of trial and channel removal
row=0;
for t = 1:length(trialind)
    
    % remove trials
    rmtrial = invvar;
    rmtrial(:,trialind(1:t)) = []; 
    
    % order variance over trials
    [~,chanind]=sort(sum(rmtrial,2),'descend');
    
    for c = 1:length(chanind)
        
        % remove chans
        rmchan = rmtrial;
        rmchan(chanind(1:c),:) = []; 
        
        % if all low var removed, calculate area left
        if sum(rmchan>1)==0
            row=row+1;
            results(row,:) = [t,c,length(rmchan(:))];
            idx{row} = {trialind(1:t),chanind(1:c)};
            continue
        end
    end
end

trial_weight = S.prep.clean.flatchan.trial_weight;
chan_weight = S.prep.clean.flatchan.chan_weights;
% use area per lost channel/trial
metric=[];
for cw = 1:length(chan_weight)
    metric(:,cw) = results(:,3)./((results(:,2)*chan_weight(cw)) + (results(:,1)*trial_weight));
end
f=figure('units','normalized','outerposition',[0 0 1 1])
for m = 1:size(metric,2)
    [~,order] = sort(metric(:,m),'descend');
    ordered_results = results(order,:);
    ordered_idx=idx(order);
    rmvar = zeros(size(invvar));
    rmvar(ordered_idx{1}{2},:) = 1;
    rmvar(:,ordered_idx{1}{1}) = 1;
    subplot(ceil(length(chan_weight)/3),ceil(length(chan_weight)/4),m);
    imagesc(rmvar)
    title(num2str(m));
    xlabel(['bad trials = ' num2str(ordered_results(1,1))]);
    ylabel(['bad chans = ' num2str(ordered_results(1,2))]);
end

if length(chan_weight)>1
    answer = str2double(inputdlg('Choose plot number','',1,{'9'}));
else
    answer = 1;
end
[~,order] = sort(metric(:,answer),'descend');
ordered_results = results(order,:);
ordered_idx=idx(order);
S.prep.clean.flatchan.rejchan = sort(ordered_idx{1}{2});
S.prep.clean.flatchan.rejtrial = sort(ordered_idx{1}{1});
close(f)
