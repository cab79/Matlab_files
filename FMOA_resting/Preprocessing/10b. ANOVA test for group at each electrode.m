
%%%%%%%%%%%%% rereferenced data %%%%%  ANOVA with all electrodes,%%%%%%%

clear all
% INPUT Args:
subjects = {'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13' 'H14', 'H15', 'H16', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'OA1', 'OA2', 'OA4', 'OA5', 'OA6', 'OA7', 'OA8', 'OA9', 'OA10', 'OA11', 'OA12', 'OA13', 'OA14', 'OA15', 'OA16', 'OA17'};
NSub = length(subjects);
Ngp1 = 16;
Ngp2 = 16;
Ngp3 = 16;
ele = [1:2 4:30 33:64]; %all electrodes
%ele = [64 28 60 56 23 17 10 30 62 29 25 24 57 61 52]; % gamma1
%ele = [64 30 62 29 25 24 57 61 28 60 23 17 10 52 56 59 63]%gamma2
%ele = [64 30 62 29 25 24 57 61 56 28 60 23 17 10 52] %gamma3
%ele = [30 62 28 60 29 25 61 57 59 10 24 64 17 23 56] %beta
%ele = [2 45 1 44 38 34 7 11 46 33 41] %alpha1
%ele = [2 45 1 44 38 7 33 41] %alpha2
%ele = [62 30 64 28 60 23 29 2 45 1 44] %theta
%ele = [64 30 62 29 28 60 23] %delta
results = zeros(1,length(ele));

%%%%%%%%%%

subject = [1:NSub 1:NSub]; %- subject ID
group = [ones(1,Ngp1) 2*ones(1,Ngp2) 3*ones(1,Ngp3)]; %-- group ID for every datapoint

%load('gamma1_closed.mat'); 
%data = gamma1_closed;
%load('gamma1_open.mat') 
%data = gamma1_open;
%load('gamma2_closed.mat')
%data = gamma2_closed;
%load('gamma2_open.mat')
%data = gamma2_open;
 %load('gamma3_open.mat') 
 %data = gamma3_open;
 %load('gamma3_closed.mat') 
 %data = gamma3_closed;
%load('beta_open.mat')
%data = beta_open;
% load('beta_closed.mat') 
 %data = beta_closed;
 %load('alpha1_open.mat') 
 %data = alpha1_open;
 load('alpha1_closed.mat')
 data = alpha1_closed;
%load('alpha2_open.mat')
%data = alpha2_open;
%load('alpha2_closed.mat')
%data = alpha2_closed;
%load('alpha3_open.mat')
%data = alpha3_open;
%load('alpha3_closed.mat')
%data = alpha3_closed;
 %load('theta_open.mat')  
 %data = theta_open;
 %load('theta_closed.mat')  
 %data = theta_closed;
%load('delta_open.mat')
%data = delta_open;
 %load('delta_closed.mat') 
%data = delta_closed;

ref_on = 1;
if ref_on==1
    ref_ele = [1:2 4:30 33:64];
    %ref_ele = [1 2 33 36 44 45 42 34 38];

    refdata = mean(data(ref_ele,:),1);
    data = data - repmat(refdata,size(data,1),1);
end

for e = 1:length(ele)

datum = data(ele(e),:);


neednan = Ngp1-Ngp3;
grp = [1 2 3];
ii1 = find(group == 1);
ii2 = find(group == 2);
ii3 = find(group == 3);
ii3 = [ii3 max(ii3)+1:max(ii3)+neednan];

datum(1,NSub+1:NSub+neednan) = NaN;
dat(:,1) = datum(ii1)';
dat(:,2) = datum(ii2)';
dat(:,3) = datum(ii3)';
[p,table,stats] = ANOVA1(dat,group,'off');
results(1,e) = p;

end

q=0.05;
%%%% FDR correction
%[pID,pN] = fdr(p,q)
%sig = find(results<=pID)
%if ~isempty(sig)
%    results(sig)
%end

%%%% no FDR correction
sig = find(results<=q)
if ~isempty(sig)
    results(sig)
end




%[elec,x,y,chans] = textread('C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs','%f %f %f %s');
%elec=elec(ele);
%ele2 = zeros(elec,1);
%chans2 = cell(elec,1);

%inter = find(results(1,:) <= pID);
%channel = chans2(inter)
%inter_pval = p(inter)
