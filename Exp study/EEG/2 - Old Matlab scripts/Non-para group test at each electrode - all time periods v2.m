clear all

subjects = {'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13' 'H14', 'H15', 'H16', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'OA1', 'OA2', 'OA4', 'OA5', 'OA6', 'OA7', 'OA8', 'OA9', 'OA10', 'OA11', 'OA12','OA13'};
NSub = length(subjects);

Ngp1 = 16;
Ngp2 = 16;
Ngp3 = 12;

ele = [1:2 4:30 33:64]; % numbers from 64 electrodes
%ele = [49 50 48 39 14 15 54 51 20]; % 9 ele around Cz
%ele = [4 14 15 20 21]; % Pz Cpz Cz FCz Fz
%ele = [14 15 20]; % Cpz Cz FCz

conditions = 6;
cond_inc = [1:6]; % conditions to include in test 


%%%%%%%%%%
fnames={
    '_1_p2_av.dat' '_2_p2_av.dat' '_3_p2_av.dat' '_4_p2_av.dat' '_5_p2_av.dat' '_6_p2_av.dat'
    '_1_n2_av.dat' '_2_n2_av.dat' '_3_n2_av.dat' '_4_n2_av.dat' '_5_n2_av.dat' '_6_n2_av.dat'
    '_1_l_av.dat' '_2_l_av.dat' '_3_l_av.dat' '_4_l_av.dat' '_5_l_av.dat' '_6_l_av.dat'
    '_1_m_av.dat' '_2_m_av.dat' '_3_m_av.dat' '_4_m_av.dat' '_5_m_av.dat' '_6_m_av.dat'
    '_1_e_av.dat' '_2_e_av.dat' '_3_e_av.dat' '_4_e_av.dat' '_5_e_av.dat' '_6_e_av.dat'
  };

peaks = size(fnames, 1);
results = zeros(length(ele),conditions,peaks);
%subject = [1:NSub 1:NSub]; %- subject ID
%time = [ones(1,NSub) 2*ones(1,NSub)]; %- timepoint
group = [ones(1,Ngp1) 2*ones(1,Ngp2) 3*ones(1,Ngp3)]; %-- group ID for every datapoint

for f = 1:peaks

    fname = fnames(f,:)';
    
for e = 1:length(ele)

data1 = zeros(1,NSub); %cond1
data2 = zeros(1,NSub); %cond2
data3 = zeros(1,NSub); %cond3
data4 = zeros(1,NSub); %cond4
data5 = zeros(1,NSub); %cond5
data6 = zeros(1,NSub); %cond6

for n = 1:NSub
sub = subjects(n);
c1 = load([char(sub) char(fname(1))]);
c2 = load([char(sub) char(fname(2))]);
c3 = load([char(sub) char(fname(3))]);
c4 = load([char(sub) char(fname(4))]);
c5 = load([char(sub) char(fname(5))]);
c6 = load([char(sub) char(fname(6))]);

c1 = c1(1,ele);
c2 = c2(1,ele);
c3 = c3(1,ele);
c4 = c4(1,ele);
c5 = c5(1,ele);
c6 = c6(1,ele);

data1(1,n) = c1(1,e);
data2(1,n) = c2(1,e);
data3(1,n) = c3(1,e);
data4(1,n) = c4(1,e);
data5(1,n) = c5(1,e);
data6(1,n) = c6(1,e);

end

% normalise the data to mean across conditions, for each subject
data4mean = reshape(cat(2,data1,data2,data3,data4,data5,data6), NSub,length(cond_inc));
for i=1:NSub
    meandata = mean(data4mean(i,:),2);
    data4mean(i,:) = data4mean(i,:) - meandata;
end
data4mean = reshape(data4mean, 1, NSub*length(cond_inc));
data1 = data4mean(1,1:NSub);
data2 = data4mean(1,NSub+1:2*NSub);
data3 = data4mean(1,2*NSub+1:3*NSub);
data4 = data4mean(1,3*NSub+1:4*NSub);
data5 = data4mean(1,4*NSub+1:5*NSub);
data6 = data4mean(1,5*NSub+1:6*NSub);

neednan = Ngp1-Ngp3;
grp = [1 2 3];
ii1 = find(group == 1);
ii2 = find(group == 2);
ii3 = find(group == 3);
ii3 = [ii3 max(ii3)+1:max(ii3)+neednan];

for c = cond_inc
data = eval(['data' num2str(c)]);
data(1,NSub+1:NSub+neednan) = NaN;
dat(:,1) = data(ii1)';
dat(:,2) = data(ii2)';
dat(:,3) = data(ii3)';
[p,table,stats] = kruskalwallis(dat,grp,'off');
results(e,c,f) = p;
end


end

end

results_pval = results(:,cond_inc,:); % '3' for t values - use '2' for n values
p = reshape(results_pval,1,length(ele)*length(cond_inc)*peaks);
q=0.05;

%%%% FDR correction
%[pID,pN] = fdr(p,q)


%% readout of significant electrodes

peak = 3;

results_pval2 = results_pval(:,:,peak);
p2 = p(1+(length(p)/peaks)*(peak-1):(length(p)/peaks)+(length(p)/peaks)*(peak-1));
p2_1 = p2(1:length(ele));
p2_2 = p2(length(ele)+1:2*length(ele));
p2_3 = p2(2*length(ele)+1:3*length(ele));
p2_4 = p2(3*length(ele)+1:4*length(ele));
p2_5 = p2(4*length(ele)+1:5*length(ele));
p2_6 = p2(5*length(ele)+1:6*length(ele));


if exist('pID') == 1
s1 = find(p2_1<=pID);
s2 = find(p2_2<=pID);
s3 = find(p2_3<=pID);
s4 = find(p2_4<=pID);
s5 = find(p2_5<=pID);
s6 = find(p2_6<=pID);
else
s1 = find(p2_1<=q);
s2 = find(p2_2<=q);
s3 = find(p2_3<=q);
s4 = find(p2_4<=q);
s5 = find(p2_5<=q);
s6 = find(p2_6<=q);
pID = q;
end

[elec,x,y,chans] = textread('C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs','%f %f %f %s');
chans = chans(ele);

cond1_chans = chans(s1)
cond1_pval = p2_1(s1)

cond2_chans = chans(s2)
cond2_pval = p2_2(s2)

cond3_chans = chans(s3)
cond3_pval = p2_3(s3)

cond4_chans = chans(s4)
cond4_pval = p2_4(s4)

cond5_chans = chans(s5)
cond5_pval = p2_5(s5)

cond6_chans = chans(s6)
cond6_pval = p2_6(s6)

