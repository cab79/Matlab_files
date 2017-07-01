clear all
close all
filepath = 'W:\Data\CRPS_Digit_Perception_exp1\';
grplist = [39 40 41 42]; 
sublist_side = {'L','R','L','R'}; %Affected vs unaffected exp1
aff_side = [2 1 1 2 1 1 1 2 1 2 1 1 1];


loadsubj
subjects = subjlists(grplist);

sublist = {};
for s = 1:length(subjects)
    for s2 = 1:length(subjects{s,1}) 
        subj = {[subjects{s,1}{s2,1} '.c1.sdnorm1'] [subjects{s,1}{s2,1} '.c5.sdnorm1']};
        
        statmode = 'trial';
        subjinfo = subj;
        condlist = {'1' '5'};
        latency = [0 0.2];
        alpha = 0.05;
        numrand = 1000; 
        ttesttail = 0;
        testgfp = 'off';
        singlesource = 'off';
        testmean = 'off';
        testlat = 'off';
        cov = [];

        stat = compamp(statmode,subjinfo,condlist,latency,cov,'alpha',alpha,'numrand',numrand,'ttesttail',ttesttail,'testgfp',testgfp,'singlesource',singlesource,'testmean',testmean,'testlat',testlat);
        pause
    end
end