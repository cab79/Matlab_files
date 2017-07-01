
clear all; clear classes
close all
filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception';
cd(filepath);
save_folder = 'C:\Data\CRPS-DP\CRPS_Digit_Perception\Results\Source_cluster';
grplist = [29]; 
%grplists = {39; 40; 41; 42};
sublist_side = {'L','R','L','R'}; %Affected vs unaffected exp1
aff_side = [2 1 1 2 1 1 1 2 1 2 1 1 1];
covariate = [];
load cov_RT_patient
%covariate = cov;
covnme = 'Acc';

regions = {%'L_postcen',60;
	'R_postcen',61;
	%'L_sup_par',62;
	%'R_sup_par',63;
	%'L_inf_par',32;
	'R_inf_par',33;
    %'peak4clus',[];
    };
Nreg = size(regions,1);

loadsubj
subjects = subjlists(grplist);

for r = 2%:Nreg

    Ns=0;
    for s = 1:length(subjects)
        for s2 = 1:length(subjects{s,1}) 
            tmp_nme = subjects{s,1}{s2,1};
            tmp_nme = strrep(tmp_nme, '.left', '_left');
            tmp_nme = strrep(tmp_nme, '.Left', '_left');
            tmp_nme = strrep(tmp_nme, '.right', '_right');
            tmp_nme = strrep(tmp_nme, '.Right', '_right');
            tmp_nme = strrep(tmp_nme, '.flip', '_flip');
            tmp_nme = strrep(tmp_nme, '.Flip', '_flip');
            tmp_nme = strrep(tmp_nme, '.aff', '_aff');
            tmp_nme = strrep(tmp_nme, '.Aff', '_aff');
            tmp_nme = strrep(tmp_nme, '.Unaff', '_unaff');
            tmp_nme = strrep(tmp_nme, '.unaff', '_unaff');
            tmp_nme = strrep(tmp_nme, '_Left', '_left');
            tmp_nme = strrep(tmp_nme, '_Right', '_right');
            tmp_nme = strrep(tmp_nme, '_Flip', '_flip');
            tmp_nme = strrep(tmp_nme, '_Aff', '_aff');
            tmp_nme = strrep(tmp_nme, '_Unaff', '_unaff');
            tmp_nme = strrep(tmp_nme, '.Exp1', '_Exp1');
            fnames{s2,s} = ['spm8_' tmp_nme '_' regions{r,1}];
            pnames{s2,s} = ['spm8_' tmp_nme];
        end
    end

    statmode = 'subj_corr';
    subjinfo = fnames;
    plotdata = pnames;
    condlist = {'Aff' 'cov'};
    condcont = [1 1]; % contrast. All '1' mean common condition values will be collapsed. '1 -1' will subtract common condition values.
    latency = [0 0.6];
    peakdef =  [1 1];
    alpha = 0.05;
    numrand = 1000; 
    ttesttail = 0;
    testgfp = 'off';
    singlesource = 'on';
    testmean = 'off';
    testlat = 'off';
    timeshift = 0.2
    %cov = covariate;

    stat = FTstats(statmode,subjinfo,condlist,condcont,latency,cov,filepath,'alpha',alpha,'numrand',numrand,'ttesttail',ttesttail,'testgfp',testgfp,...
    'singlesource',singlesource,'testmean',testmean,'testlat',testlat,'timeshift',timeshift,'peakdef',peakdef);

    save_name = [regions{r,1} '_grp' num2str(grplist) '_' covnme];
    save(fullfile(save_folder,save_name),'stat');
    h = get(0,'children');
    for i=1:length(h)
        saveas(h(i), fullfile(save_folder,[save_name '_' get(h(i),'Name')]), 'fig');
    end
end