clear all

subjects = {'H1','F9', 'OA8'};

fnames={
    '_1_avg_ca.mat';
    '_2_avg_ca.mat';    
    '_3_avg_ca.mat';
    '_4_avg_ca.mat';
    '_5_avg_ca.mat';
    '_6_avg_ca.mat'; 
  };

for n = 2
    subject = subjects(n);
    subject = char(subject);
for x=1:length(fnames);
    fname = fnames(x);
    fname = char(fname);
    sub_fname= [subject fname];
    load(sub_fname); 
    
    avg(36,:) = mean(avg([35 37],:),1);
    
    eval(['save ' sub_fname ' avg']);
end
    
end