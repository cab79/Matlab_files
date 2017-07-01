function [freqs_nme,freqs_limits,freq_idx] = find_freq_limits(gavg,no_freqs)
% gavg is a call array of grand averages, size (num_avgs,1)

nav = size(gavg,1);

% cross-correlations
cormat=[];
for s = 1:nav 
    for e = 1:size(gavg{s,1},1)
        A = squeeze(gavg{s,1}(e,:,:));
        cormat(e,:,:,s) = corrcoef(A(~any(isnan(A),2),:));
    end
end
cormatm = mean(cormat,4);
corav = squeeze(mean(cormatm,1));

% Freq bin identification
for threshcor=0.8:0.001:1;
    corth=corav.*(corav>threshcor);
    bins = conncomp(graph(corth));
    if length(unique(bins))==no_freqs
        break
    end
end

freq_idx={};
for f = 1:no_freqs
    freq_idx{f}=find(bins==f);
end
%for cc=1:size(corth,2)
%    if any(corth(:,cc))==0
%        freq_idx{length(freq_idx)+1}=cc;
%    elseif any(corth(:,cc))
%        tempi=find(corth(:,cc)==1);
%        if length(freq_idx)>1
%            [te fr] = intersect(tempi,1:max(freq_idx{cc-1}));
%            tempi(fr)=[];
%        end
%        freq_idx{length(freq_idx)+1}=tempi;
%        if any(tempi==size(corth,1))
%            break
%        end
%    end
%end

freqs_nme = cell(no_freqs,1);
freqs_limits = cell(no_freqs,1);
for f = 1:no_freqs
    fr = freqs(freq_idx{f});
    C1 = strsplit(num2str(fr(1)),'.');
    if length(fr)>1
        C2 = strsplit(num2str(fr(end)),'.');
        freqs_nme{f} = ['f' C1{1} '_f' C2{1}];
    else
        freqs_nme{f} = ['f' C1{1}];
    end
    freqs_limits{f} = [fr(1) fr(end)];
end