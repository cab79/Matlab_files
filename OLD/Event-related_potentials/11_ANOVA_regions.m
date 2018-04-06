clear all

chan_locs = 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan_MR62.locs';
[tmpeloc labels Th Rd indices] = readlocs(chan_locs,'filetype','loc');
labels = upper(labels);
p1 = find(ismember(labels,{'AF7','F3','F5','F7'}));
p2 = find(ismember(labels,{'FZ','F1','F3','AFZ','AF1','AF3','FPZ','FP1','FP3'}));
p3 = find(ismember(labels,{'AF8','F4','F6','F8'}));
p4 = find(ismember(labels,{'FT7','FC5','FC3','T7','C5','C3','TP7','CP5','CP3'}));
p5 = find(ismember(labels,{'FC1','FC2','CZ','C1','C2','CPZ','CP1','CP2'}));
p6 = find(ismember(labels,{'FT8','FC6','FC4','T8','C6','C4','TP8','CP6','CP4'}));
p7 = find(ismember(labels,{'P7','P3','P5','PO7'}));
p8 = find(ismember(labels,{'PZ','P1','P2','POZ','PO3','PO4','OZ','O1','O2'}));
p9 = find(ismember(labels,{'P8','P4','P6','PO8'}));
region_labels = {'anterior left';'anterior mid'; 'anterior right';'central left';'central mid'; 'central right';'posterior left';'posterior mid'; 'posterior right';};
files = {'delta_9reg.mat';
    'theta_9reg.mat';
    'alpha_9reg.mat';
    'beta_9reg.mat'; 
    'gamma1_9reg.mat';
    'gamma2_9reg.mat'
}; 

Ptable = cell(1,1);
Var = {'region','scan','time','region*scan','region*time','scan*time','region*scan*time'};
Ptable(2:8,1) = Var;

for f = 3%1:length(files)
    
load(files{f});
[pth,fname,ext] = fileparts(files{f});
Ptable{1,f+1} = fname;

% Change data dimensions to [subject,region,scan,time]
ndata = permute(ndata,[2 1 3 4]);
sd = size(ndata);
% reshape such that time is the largest factor
data = reshape(ndata,1,sd(1)*sd(2)*sd(3)*sd(4));
% create design matrix 
subjects = repmat([1:sd(1)],1,sd(2)*sd(3)*sd(4));
regions = repmat(sort(repmat([1:sd(2)],1,sd(1)), 'ascend'),1,sd(3)*sd(4));
scans = repmat(sort(repmat([1:sd(3)],1,sd(1)*sd(2)), 'ascend'),1,sd(4));
time = sort(repmat([1:sd(4)],1,sd(1)*sd(2)*sd(3)), 'ascend');
design = [data' regions' scans' time' subjects'];

%[p,table,stats,terms] = anovan(design(:,1),{subjects' regions' scans' time'},'model','full');

[Pv] = RMAOV33(design,0.05);
Ptable(2:8,f+1) = Pv';

end

%ANOVA for each time, one frequency only
P1table = cell(1,1);
Var = {'scan','region','scan*region'};
P1table(2:4,1) = Var;
P1table(1,2:4) = {'Pre', 'Pain', 'Post'};
for t = 1:3
design1 = design(design(:,4)==t,:);
[Pv1] = RMAOV2_mod(design1(:,[1 2 3 5]),0.05);
P1table(2:4,t+1) = num2cell(Pv1');
end

%ANOVA for each region, one frequency only
P2table = cell(1,1);
Var = {'scan','time','scan*time'};
P2table(2:4,1) = Var;
P2table(1,2:10) = region_labels;
for r = 1:9
design2 = design(design(:,2)==r,:);
[Pv2] = RMAOV2_mod(design2(:,[1 3 4 5]),0.05);
P2table(2:4,r+1) = num2cell(Pv2');
end

%ANOVA for each scan, one frequency only
P3table = cell(1,1);
Var = {'region','time','region*time'};
P3table(2:4,1) = Var;
P3table(1,2:3) = {'pain','no pain'};
for s = 1:2
design3 = design(design(:,3)==s,:);
[Pv3] = RMAOV2_mod(design3(:,[1 2 4 5]),0.05);
P3table(2:4,s+1) = num2cell(Pv3');
end

%ANOVA of time for each scan and region, one frequency only
P4table = cell(1,1);
Var = {'pain','no pain'};
P4table(2:3,1) = Var';
P4table(1,2:10) = region_labels;
%P4table(1,2:3) = {'pain','no pain'};
for s = 1:2
    for r = 1:9
        design4 = design(design(:,3)==s,:);
        design4 = design4(design4(:,2)==r,:);
        Pv4 = RMAOV1_mod(design4(:,[1 4 5]),0.05);
        P4table(1+s,1+r) = num2cell(Pv4');
    end
end

%ANOVA of regions for each scan and time, one frequency only
P5table = cell(1,1);
Var = {'pain','no pain'};
P5table(2:3,1) = Var';
P5table(1,2:4) = {'Pre', 'Pain', 'Post'};
for s = 1:2
    for t = 1:3
        design5 = design(design(:,3)==s,:);
        design5 = design5(design5(:,4)==t,:);
        Pv5 = RMAOV1_mod(design5(:,[1 2 5]),0.05);
        P5table(1+s,1+t) = num2cell(Pv5');
    end
end

% plot scan x time for each region separately
for r = 1:9;
figure(r)
X = design(design(:,2)==r,1);
S = design(design(:,2)==r,3);
T = design(design(:,2)==r,4);
ST = {S T};
labels = {'Pre' 'Pain' 'Post' 'Pre' 'Pain' 'Post'};
name = char(region_labels(r));
boxplot(X,ST,'plotstyle','traditional','boxstyle','outline','colorgroup',ST,'color','rrrbbb','factorseparator',[1],'factorgap',[20 1],'jitter',0,'labels',labels,'notch','off','outliersize',3,'symbol','o'), title(name);
end

% plot region x time for both scans
for s = 1:2
design3 = design(design(:,3)==s,:);
figure(s)
colorset = varycolor(3);
colorset = repmat(colorset,9,1);
X = design3(:,1);
R = design3(:,2);
T = design3(:,4);
RT = {R T};
labels = {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' ''};
name = '';
boxplot(X,RT,'plotstyle','traditional','boxstyle','outline','colorgroup',RT,'color',colorset,'factorseparator',[1],'factorgap',[20 1],'jitter',0,'labels',labels,'notch','off','outliersize',3,'symbol','o'), title(name);
end

%reshape for SPSS
data_SPSS = reshape(data,max(design(:,5)),max(design(:,3))*max(design(:,4))*max(design(:,2)));
