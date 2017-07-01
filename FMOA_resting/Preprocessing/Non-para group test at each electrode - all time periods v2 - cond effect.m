clear all

subjects = {'M1';'M2';'M4';'M6';'M7';'M10';'M12';'M14';'M15';'M16';'M17';'M29';'M30';'M32';'M33';'M9';'M19';'M20';'M21';'M22';'M24';'M25';'M26';'M35';'M36';'M37';'M38';'M40'};
NSub = length(subjects);

Ngp1 = 15;
Ngp2 = 13;

%ele = [1:2 4:30 33:64]; % numbers from 64 electrodes
ele = [49 50 48 39 14 15 54 51 20]; % 9 ele around Cz
%ele = [4 14 15 20 21]; % Pz Cpz Cz FCz Fz
%ele = [14 15 20]; % Cpz Cz FCz

conditions = 1;
cond_inc = [1]; % conditions to include in test 


%%%%%%%%%%
fnames={
    '_5_p2_av.dat' '_6_p2_av.dat' 
    '_5_l_av.dat' '_6_l_av.dat' 
  };

peaks = size(fnames, 1);
results = zeros(length(ele),conditions,peaks);
%subject = [1:NSub 1:NSub]; %- subject ID
%time = [ones(1,NSub) 2*ones(1,NSub)]; %- timepoint
group = [ones(1,Ngp1) 2*ones(1,Ngp1)]; %-- group ID for every datapoint

for f = 1:peaks

    fname = fnames(f,:)';
    
for e = 1:length(ele)

data1 = zeros(1,NSub); %cond1


for n = 1:NSub
sub = subjects(n);
c1 = load([char(sub) char(fname(1))]);
c2 = load([char(sub) char(fname(2))]);


c1 = c1(1,ele);
c2 = c2(1,ele);


data1(1,n) = c1(1,e) - c2(1,e);


end

% normalise the data to percent of mean across conditions, for each subject
%data4mean = reshape(cat(2,data1,data2), NSub,2);
%for i=1:NSub
%    meandata = mean(data4mean(i,:),2);
%    %data4mean(i,:) = ((data4mean(i,:)-meandata)/meandata)*100;
%    data4mean(i,:) = data4mean(i,:) - meandata;
%end
%data4mean = reshape(data4mean, 1, NSub*2);
%data1 = data4mean(1,1:NSub);
%data2 = data4mean(1,NSub+1:2*NSub);

neednan = Ngp1-Ngp2;
grp = [1 2];
ii1 = find(group == 1);
ii2 = find(group == 2);

data1(1,Ngp1+Ngp2+1:Ngp1+Ngp2+neednan) = NaN;
data(:,1) = data1(ii1)';
data(:,2) = data1(ii2)';
[p,table,stats] = kruskalwallis(data,grp,'off');
results(e,1,f) = p;

end

end

results_pval = results(:,cond_inc,:); % '3' for t values - use '2' for n values
p = reshape(results_pval,1,e*length(cond_inc)*peaks);
q=0.05;

%%%% FDR correction
%[pID,pN] = fdr(p,q)


%% readout of significant electrodes

peak = 2;

results_pval2 = results_pval(:,:,peak);
p2 = p(1+(length(p)/peaks)*(peak-1):(length(p)/peaks)+(length(p)/peaks)*(peak-1));

if exist('pID') == 1
sig = find(p2<=pID);
else
sig = find(p2<=q);
pID = q;
end

s1 = find(sig<=length(ele));

[elec,x,y,chans] = textread('C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs','%f %f %f %s');
chans = chans(ele);
ele2 = [ele];
chans2 = chans;

cond1_chans = chans2(sig(s1))
cond1_pval = p2(sig(s1))


