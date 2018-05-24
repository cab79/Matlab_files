function S = flat_channel_reject(S,EEG)
varthresh = S.prep.clean.flatchan.varthresh;
results = [];
idx={};

invvar = 1./squeeze(std(EEG.data,[],2));
%imagesc(invvar)

% order variance over trials
[~,trialind]=sort(sum(invvar,1),'descend');
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
        
        % for every combination of trial and channel removal, 

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

answer = str2double(inputdlg('Choose plot number'))
[~,order] = sort(metric(:,answer),'descend');
ordered_results = results(order,:);
ordered_idx=idx(order);
S.prep.clean.flatchan.rejchan = sort(ordered_idx{1}{2});
S.prep.clean.flatchan.rejtrial = sort(ordered_idx{1}{1});
close(f)
