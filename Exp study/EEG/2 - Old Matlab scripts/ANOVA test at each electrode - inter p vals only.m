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

ele = [1:2 4:30 33:64]; % numbers from 64 electrodes
%ele = [49 50 48 39 14 15 54 51 20]; % 9 ele around Cz
%ele = [49 50 48 39 14 15 54 51 20 4 35 37]; % 12 ele around Cz
%ele = [4 14 15 20 21]; % Pz Cpz Cz FCz Fz
%ele = [14 15 20]; % Cpz Cz FCz
%ele = [4 14 15 20 36]; % 5 midline ele
%ele = [49 50 48 39 14 15 4 35 37]; % 9 ele around CPz
%elec = {
 %   [34 36 38];
  %  [14 39 48];
 %   [15 49 50];
 %   [15 49 50];
   % [4 35 37]
  %  };

conditions = 1;
cond_inc = [1]; % conditions to include in test 


%%%%%%%%%%
fnames={
    '_2_e_av.dat' '_4_e_av.dat'
   % '_2_l_av.dat' '_4_l_av.dat'
  %  '_2_n2_av.dat' '_4_n2_av.dat'
  % '_2_p2_av.dat' '_4_p2_av.dat'
   % '_2_po_av.dat' '_4_po_av.dat'
  };

peaks = size(fnames, 1);
results = zeros(3,3,1,conditions,peaks);
subject = [1:NSub 1:NSub]; %- subject ID
time = [ones(1,NSub) 2*ones(1,NSub)]; %- timepoint
group = [ones(1,Ngp1) 2*ones(1,Ngp2) ones(1,Ngp1) 2*ones(1,Ngp2)]; %-- group ID for every datapoint

for f = 1:peaks

    fname = fnames(f,:)';
    if exist('elec') ==1
    ele = elec(f);
    ele = ele{:};
    end
    
for e = 1:length(ele)

data1 = zeros(1,NSub*2); %cond1
%data2 = zeros(1,NSub*2); %cond2
%data3 = zeros(1,NSub*2); %cond2-1

for n = 1:NSub
sub = subjects(n);
c1 = load([char(sub) char(fname(1))]);
c2 = load([char(sub) char(fname(2))]);
%c3 = load([char(sub) char(fname(3))]);
%c4 = load([char(sub) char(fname(4))]);
%c5 = load([char(sub) char(fname(5))]);
%c6 = load([char(sub) char(fname(6))]);

c1 = c1(1,ele);
c2 = c2(1,ele);
%c3 = c3(1,ele);
%c4 = c4(1,ele);
%c5 = c5(1,ele);
%c6 = c6(1,ele);

data1(1,n) = c1(1,e);
%data2(1,n) = c2(1,e);
%data3(1,n) = c5(1,e);
data1(1,n+NSub) = c2(1,e);
%data2(1,n+NSub) = c4(1,e);
%data3(1,n+NSub) = c6(1,e);

end

% normalise the data to percent of mean across conditions, for each subject
data4mean = reshape(data1, NSub,2);
for i=1:NSub
    meandata = mean(data4mean(i,:),2);
    %data4mean(i,:) = ((data4mean(i,:)-meandata)/meandata)*100;
    data4mean(i,:) = data4mean(i,:) - meandata;
end
data4mean = reshape(data4mean, 1, NSub*2);
data1 = data4mean(1,1:2*NSub);
%data2 = data4mean(1,2*NSub+1:4*NSub);

data = data1;
tab = F1LDF1(data,group,subject,time);
results(:,:,e,1,f) = tab;

%data = data2;
%tab = F1LDF1(data,group,subject,time);
%results(:,:,e,2,f) = tab;

%data = data3;
%tab = F1LDF1(data,group,subject,time);
%results(:,:,e,3,f) = tab;


end

end

results_pval = results(:,3,:,cond_inc,:); % '3' for t values - use '2' for n values
results_pval = results(3,3,:,cond_inc,:); % interaction p values only
p = reshape(results_pval,1,e*length(cond_inc)*peaks);
q=0.05;

%%%% FDR correction
%[pID,pN] = fdr(p,q)


%% readout of significant electrodes

peak = 1;
if exist('elec') ==1
ele = elec(peak);
ele = ele{:};
end

results_pval2 = results_pval(:,:,:,:,peak);
p2 = p(1+(length(p)/peaks)*(peak-1):(length(p)/peaks)+(length(p)/peaks)*(peak-1));

if exist('pID') == 1
sig = find(p2<=pID);
else
sig = find(p2<=q);
pID = q;
end

s1 = find(sig<=length(ele));
s2 = find(sig<=length(ele)*length(cond_inc));
if length(cond_inc) > 2
s3 = find(sig<=length(ele)*length(cond_inc));
s3(s2) = [];
end
s2(s1) = [];

[elec,x,y,chans] = textread('C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs','%f %f %f %s');
chans = chans(ele);
ele2 = ele
chans2 = chans

cond1_chans = chans2(sig(s1))
c1p = ismember(ele,ele2(sig(s1)));
cond1_pval = results_pval2(:,:,c1p,1)
interp1 = find(cond1_pval(:,:,:) <= pID)
inter1 = cond1_chans(interp1)
