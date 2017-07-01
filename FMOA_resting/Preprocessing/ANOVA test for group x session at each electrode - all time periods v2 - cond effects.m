clear all

% INPUT Args:
subjects = {'M1';'M2';'M4';'M6';'M7';'M10';'M12';'M14';'M15';'M16';'M17';'M29';'M30';'M32';'M33';'M9';'M19';'M20';'M21';'M22';'M24';'M25';'M26';'M35';'M36';'M37';'M38';'M40'};
% reduced n
%subjects = {'M1';'M2';'M6';'M7';'M10';'M14';'M15';'M16';'M17';'M29';'M32';'M33';'M20';'M21';'M22';'M24';'M25';'M26';'M35';'M36';'M37';'M38';'M40'};
NSub = length(subjects);

Ngp1 = 15;
Ngp2 = 13;

%reduced n
%Ngp1 = 12;
%Ngp2 = 11;

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
results = zeros(3,3,length(ele),conditions,peaks);
subject = [1:NSub 1:NSub]; %- subject ID
time = [ones(1,NSub) 2*ones(1,NSub)]; %- timepoint
group = [ones(1,Ngp1) 2*ones(1,Ngp2) ones(1,Ngp1) 2*ones(1,Ngp2)]; %-- group ID for every datapoint

for f = 1:peaks

    fname = fnames(f,:)';
    
for e = 1:length(ele)

data3 = zeros(1,NSub*2); %cond2-1

for n = 1:NSub
sub = subjects(n);
c5 = load([char(sub) char(fname(1))]);
c6 = load([char(sub) char(fname(2))]);

c5 = c5(1,ele);
c6 = c6(1,ele);

data3(1,n) = c5(1,e);
data3(1,n+NSub) = c6(1,e);

end

% normalise the data to percent of mean across conditions, for each subject
%data4mean = reshape(cat(2,data1,data2), NSub,4);
%for i=1:NSub
%    meandata = mean(data4mean(i,:),2);
%    %data4mean(i,:) = ((data4mean(i,:)-meandata)/meandata)*100;
%    data4mean(i,:) = data4mean(i,:) - meandata;
%end
%data4mean = reshape(data4mean, 1, NSub*4);
%data1 = data4mean(1,1:2*NSub);
%data2 = data4mean(1,2*NSub+1:4*NSub);

data = data3;
tab = F1LDF1(data,group,subject,time);
results(:,:,e,1,f) = tab;

end
end

results_pval = results(:,3,:,cond_inc,:); % '3' for t values - use '2' for n values
p = reshape(results_pval,1,e*length(cond_inc)*3*peaks);
q=0.05;

%%%% FDR correction
%[pID,pN] = fdr(p,q)


%% readout of significant electrodes

peak = 2;

results_pval2 = results_pval(:,:,:,:,peak);
p2 = p(1+(length(p)/peaks)*(peak-1):(length(p)/peaks)+(length(p)/peaks)*(peak-1));

if exist('pID') == 1
sig = find(p2<=pID);
else
sig = find(p2<=q);
pID = q;
end

s3 = find(sig<=3*length(ele));

[elec,x,y,chans] = textread('C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs','%f %f %f %s');
chans = chans(ele);
ele2 = zeros(length(ele)*length(cond_inc)*3,1);
chans2 = cell(length(ele)*length(cond_inc)*3,1);
for c = 1:length(ele)
chans2(1+3*(c-1):3+3*(c-1)) = chans(c);
chans2(length(ele)*3+1+3*(c-1):length(ele)*3+3+3*(c-1)) = chans(c);
chans2(2*length(ele)*3+1+3*(c-1):2*length(ele)*3+3+3*(c-1)) = chans(c);
ele2(1+3*(c-1):3+3*(c-1)) = ele(c);
ele2(length(ele)*3+1+3*(c-1):length(ele)*3+3+3*(c-1)) = ele(c);
ele2(2*length(ele)*3+1+3*(c-1):2*length(ele)*3+3+3*(c-1)) = ele(c);
end

cond3_chans = chans2(sig(s3))
c3p = ismember(ele,ele2(sig(s3)));
cond3_pval = results_pval2(:,:,c3p)
interp3 = find(cond3_pval(3,:,:) <= pID)
inter3 = cond3_chans(interp3)
