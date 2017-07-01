
clear all; clear classes
close all
filepath = 'W:\Data\CRPS_Digit_Perception_exp1\';
cd(filepath);
save_folder = 'W:\Data\CRPS_Digit_Perception_exp1\Results\Fun_conn';
%filepath = '/scratch/cb802/Data/CRPS_Digit_Perception_exp1/';
grplist = [39 40 41 42]; 
%grplists = {39; 40; 41; 42};
sublist_side = {'L','R','L','R'}; %Affected vs unaffected exp1
aff_side = [2 1 1 2 1 1 1 2 1 2 1 1 1];

%regions = {
%    'Rpar',[];
%    'Lpar',[];
%    'Rfron',[];
%    'Lfron',[];
%    'R_postcen',[];
 %   'L_postcen',[];
 %   'S1',[];
 %   };
 
 regions = {
    'MSP_Rinfpar',[];
    'MSP_Rsuppar',[];
    'MSP_Lsuppar',[];
    };

Nreg = size(regions,1);

%corrmat = [
%    1 0 1 0 0 0 0;
%    1 0 0 1 0 0 0;
%    0 1 1 0 0 0 0;
%    0 1 0 1 0 0 0;
%    1 0 0 0 1 0 0;
%    1 0 0 0 0 0 1;
%    0 1 0 0 0 1 0;
%    0 1 0 0 0 0 1;
%    0 0 1 0 1 0 0;
%    0 0 1 0 0 0 1;
%    0 0 0 1 0 1 0;
%    0 0 0 1 0 0 1;
%    ];

corrmat = [
    1 1 0;
    1 0 1;
    0 1 1;
    ];

loadsubj
subjects = subjlists(grplist);

for c = 1:size(corrmat,1)
    ri = find(corrmat(c,:));
    corrdat{1,c+1} = [regions{ri(1),1} '_' regions{ri(2),1}]
end

Ns=0;
ind=0;
for s = 1:length(subjects)
    for s2 = 1:length(subjects{s,1}) 
        ind = ind+1;
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
        corrdat{ind+1,1} = tmp_nme;
        for r = 1:Nreg
            fname = ['spm8_' tmp_nme '_' regions{r,1}];
            D = spm_eeg_load(fullfile(filepath,fname));
            F = spm2fieldtrip(D);
            dat = [];
            for f = 1:length(F.trial)
                dat = [dat F.trial{1,f}];
            end
            ROIdat(:,r) = dat';
        end
        
        for c = 1:size(corrmat,1)
            ri = find(corrmat(c,:));
            corrdat{ind+1,c+1} = corr(ROIdat(:,ri(1)),ROIdat(:,ri(2)));
        end
    end
end

save(fullfile(save_folder,'MSP_par.mat'), 'corrdat');
