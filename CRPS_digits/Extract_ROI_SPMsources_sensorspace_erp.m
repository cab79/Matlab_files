clear all
close all
dname = pwd;

cd(dname);
grplist = [39 40 41 42];  %Affected vs unaffected exp1
%grplist = [33 34 31 32]; %Affected vs %unaffected exp2

no_cond = 5; % no of conditions per data file (arm)

Rdir = 'C:\Data\CRPS-DP\CRPS_Digit_Perception_exp1\correcttrials\SPM image files\Sensorspace_masks';
%Rdata = 'W:\Brain_atlas_Hammers2003\Hammers_mith_atlas_n30r83_SPM5.img';

regions = {'Exp1_268.nii';
    };
Nreg = size(regions,1);

use_flipped = 1;
use_flip_reorient = 0;
suffix='';

%filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception\SPM image files\eleposflip\Sensorspace_images';
filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception_exp1\alltrials\SPM image files\eleposflip\Sensorspace_images';
run('C:\Data\Matlab\Matlab_files\CRPS_digits\loadsubj.m');

subjects = subjlists(grplist);

results = cell(1,length(subjects));

for s = 1:length(subjects)
    results{1,s} = cell(1,length(subjects{s,1}));
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
        
        if use_flipped==1
            if strfind(tmp_nme, 'left')
                tmp_nme = ['maspm8_' tmp_nme suffix];
                trials = [1:5];
            elseif strfind(tmp_nme, 'right')
                tmp_nme = ['maspm8_flip_' tmp_nme suffix];
                trials = [6:10];
            end
        else
            if strfind(tmp_nme, 'left')
                trials = [1:5];
            elseif strfind(tmp_nme, 'right')
                trials = [6:10];
            end
            tmp_nme = ['maspm8_' tmp_nme suffix];
        end
        
        fnames = [];
        for i = 1:no_cond
            ind = (s2-1)*no_cond + i;
            if ~isempty(strfind(tmp_nme, 'right')) && use_flip_reorient
                nme = dir(fullfile(filepath,tmp_nme,['type_' num2str(trials(i))],'strial*_flip_reorient.img'));
            else
                nme = dir(fullfile(filepath,tmp_nme,['type_' num2str(trials(i))],'strial*.img'));
            end
            for n = 1:length(nme)
                if (~isempty(strfind(tmp_nme, 'left')) && ~isempty(strfind(nme(n).name,'flip'))) || ~isempty(strfind(nme(n).name, 'Exp1')); continue; end;
                fnames{length(fnames)+1,1} = fullfile(filepath,tmp_nme,['type_' num2str(trials(i))],nme(n).name);
            end
        end
        
        for r = 1:Nreg
            for f = 1:length(fnames)
                Pdata = fnames{f,1};
                Rdata = fullfile(Rdir,regions{r});
                input{1,1} = Pdata;
                input{2,1} = Rdata;
                expres = 'i1.*i2';
                [pth nm ext] = fileparts(fnames{f,1});
                Pfname_out = fullfile(pth, [nm '_' regions{r}]);
                if ~exist(Pfname_out,'file') Output = spm_imcalc_ui(input,Pfname_out,expres);
                end
            end
        end

        
        results{1,s}{1,s2} = cell(length(fnames)+1,Nreg+1);
        for f = 1:length(fnames)
            ns = strfind(fnames{f,1},'type');
            results{1,s}{1,s2}{f+1,1} = str2num(strrep(fnames{f,1}(ns+5:ns+6),'\',''));
        end
        results{1,s}{1,s2}(1,2:end) = regions';
        results{1,s}{1,s2}(2:end,2:end) = num2cell(NaN(size(results{1,s}{1,s2}(2:end,2:end))));

        Rnii = load_nii(Rdata);
        for r = 1:Nreg
            region = regions{r};
            Rsize = length(find(Rnii.img==1)); 
            for f = 1:length(fnames)
                [pth nm ext] = fileparts(fnames{f,1});
                Pdata = fullfile(pth, [nm '_' regions{r}]);
                nii = load_nii(Pdata);
                nii_all = nii.img(nii.img~=0);
                nii_mean = sum(nii_all)/Rsize;
                results{1,s}{1,s2}{f+1,(1+r)} = nii_mean;
            end
        end
        [sorted si] = sort([results{1,s}{1,s2}{2:end,1}]);
        for col = 1:length(results{1,s}{1,s2}(1,:))
            results{1,s}{1,s2}(2:end,col) = results{1,s}{1,s2}(si+1,col);
        end
        
    end
end



save clustersERP.mat results;

results2 = [];
for s = 1:length(subjects)
    for s2 = 1:length(subjects{s,1}) 
        results2(size(results2,1)+1:size(results2,1)+5,1) = [results{1,s}{1,s2}{2:6,2}];
    end
end
results2 = reshape(results2,size(results2,1)/s,s);

% generate correlation coefficients
%slopes = [];
%for s = 1:length(subjects)
%    for s2 = 1:length(subjects{s,1}) 
%        for r = 1:Nreg
%            rho = corr(single([results{1,s}{1,s2}{2:end,1}]'), single([results{1,s}{1,s2}{2:end,r+1}]'),'type','Spearman');
%            slopes(r,s,s2) = rho;
%        end
%    end
%end

%close all
%for r = 1:Nreg
%    figure
%    X = 1:4;
%    Y = squeeze(mean(slopes(r,:,:),3));
%    E = std(squeeze(slopes(r,:,:))')/sqrt(size(slopes(r,:,:),3));
%    errorbar(X,Y,E);
%    p = anova1(squeeze(slopes(r,:,:))');
%end
