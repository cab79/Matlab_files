clear all
close all

grplist = [39 40 41 42];  %Affected vs unaffected exp1
%grplist = [35 36 37 38]; 
%grplist = [1 2 29 30]; %sublist_side = {'L','R','L','R'}; %Affected vs %unaffected exp2

no_cond = 5; % no of conditions per data file (arm)

Rdir = 'C:\Data\CRPS-DP\CRPS_Digit_Perception_exp1\correcttrials\SPM image files\Sensorspace_masks';
%Rdata = 'W:\Brain_atlas_Hammers2003\Hammers_mith_atlas_n30r83_SPM5.img';

regions = {
    'Exp1_correct_132.nii';
    'Exp1_268.nii';
    };
use_flipped=1;
Nreg = size(regions,1);

%if isunix
%    filepath = '/scratch/cb802/Data/CRPS_Digit_Perception/SPM image files/Sensorspace_images';
%    run('/scratch/cb802/Matlab_files/CRPS_digits/loadsubj.m');
%else
    filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception_exp1\alltrials';
   run('C:\Data\Matlab\Matlab_files\CRPS_digits\loadsubj.m');
%end

subjects = subjlists(grplist);

results = cell(1,length(subjects));
results_corr = [];

for s = 1:length(subjects)
    results{1,s} = cell(1,length(subjects{s,1}));
    for s2 = 1:length(subjects{s,1}) 
        tmp_nme = subjects{s,1}{s2,1}
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
                tmp_nme = ['aspm8_' tmp_nme];
                trials = [1:5];
            elseif strfind(tmp_nme, 'right')
                tmp_nme = ['aspm8_flip_' tmp_nme];
                trials = [6:10];
            end
        else
            if strfind(tmp_nme, 'left')
                trials = [1:5];
            elseif strfind(tmp_nme, 'right')
                trials = [6:10];
            end
            tmp_nme = ['aspm8_' tmp_nme];
        end
        
        % list all filenames for all trials and all conditions
        fnames = [];
        for i = 1:no_cond
            %ind = (s2-1)*no_cond + i;
            nme = dir(fullfile(filepath,tmp_nme,['type_' num2str(trials(i))],'strial*.img'));
            for n = 1:length(nme)
                fnames{length(fnames)+1,1} = fullfile(filepath,tmp_nme,['type_' num2str(trials(i))],nme(n).name);
            end
        end

        % create results output matrix
        results{1,s}{1,s2} = cell(length(fnames)+1,Nreg+1);
        for f = 1:length(fnames)
            results{1,s}{1,s2}{f+1,1} = fnames{f,1}(end-21:end-4);
        end
        results{1,s}{1,s2}(1,2:end) = regions';
        results{1,s}{1,s2}(2:end,2:end) = num2cell(NaN(size(results{1,s}{1,s2}(2:end,2:end))));
        

        % 
        %for r = 1:Nreg
        %    for f = 1:length(fnames)
        %        Pdata = fnames{f,1};
        %        Rdata = fullfile(Rdir,regions{r});
        %        input{1,1} = Pdata;
        %        input{2,1} = Rdata;
        %        expres = 'i1.*i2';
        %        [pth nm ext] = fileparts(fnames{f,1});
        %        Pfname_out = fullfile(pth, [nm '_' regions{r}]);
        %        if ~exist(Pfname_out,'file') Output = spm_imcalc_ui(input,Pfname_out,expres);
        %        end
        %    end
        %end
        
        
        for r = 1:Nreg
            Rdata = fullfile(Rdir,regions{r});
            Rnii = load_nii(Rdata);
            Rsize = length(find(Rnii.img==1)); 
            for f = 1:length(fnames)
                %[pth nm ext] = fileparts(fnames{f,1});
                %Pdata = fullfile(pth, [nm '_' regions{r}]);
                nii = load_nii(fnames{f,1});
                nii.img = permute(nii.img,[2 1 3]); % CHECK REG FILE AND DATA ARE THE SAME ORIENTATION; IF NOT CHANGE DATA ORIENTATION HERE
                roi = nii.img.*single(Rnii.img);
                roi_mean = nansum(roi(:))/Rsize;
                results{1,s}{1,s2}{f+1,(1+r)} = roi_mean;
            end
        end
        %[sorted si] = sort([results{1,s}{1,s2}{2:end,1}]);
        %for col = 1:length(results{1,s}{1,s2}(1,:))
        %    results{1,s}{1,s2}(2:end,col) = results{1,s}{1,s2}(si+1,col);
        %end
        
        % cross-correlation of ROIs
        corrmat = corr(cell2mat(results{1,s}{1,s2}(2:end,2:end)),'type','Spearman');
        
        % the following lines are specific to the case of only 2 ROIs being
        % correlated with each other as the output.
        results_corr(s,s2) = corrmat(1,2);
        
        
    end
end

save('roi_results.mat','results','results_corr');

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
