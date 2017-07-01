clear all
close all
dname = 'W:\Data\CRPS_Digit_Perception_exp1\SPM image files\LORflip-ind';

cd(dname);
grplist = [39 40 41 42]; sublist_side = {'L','R','L','R'}; %Affected vs unaffected exp1
%grplist = [33 34 31 32]; sublist_side = {'L','R','L','R'}; %Affected vs %unaffected exp2
peaks_nme = {'5'};
cond_nme = {'1','2','3','4','5'};
hand_nme = {'L','R'};
no_cond = length(cond_nme); % no of conditions per data file (arm)
no_peaks = length(peaks_nme);


Rdata = 'W:\Data\CRPS_Digit_Perception_exp1\Results\SPM source stats\Regression - indCond\flipped\LORind\Acc\subject effect\Ranked\LOR\130ms';
%Rdata = 'W:\Brain_atlas_Hammers2003\Hammers_mith_atlas_n30r83_SPM5.img';

regions = {'all_clusters.img';
    };
Nreg = size(regions,1);

loadsubj
subjects = subjlists(grplist);

Ns=0;
for s = 1:length(subjects)
    for s2 = 1:length(subjects{s,1}) 
        Ns=Ns+1;
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
        if strfind(tmp_nme,'right')
            fnames{Ns,1} = ['maspm8_flip_' tmp_nme];
        else
            fnames{Ns,1} = ['maspm8_' tmp_nme];
        end
    end
end

for r = 1:Nreg
    region = fullfile(Rdata,regions{r,1});
    [pth region_nme ext] = fileparts(region);
    for f = 1:length(fnames)
        for p = 1:no_peaks
            for c = 1:no_cond
                Pfname = fullfile(dname,[fnames{f,1} '_' peaks_nme{1,p} '_' cond_nme{1,c}]);
                Pdata = [Pfname '.nii'];
                input{1,1} = Pdata;
                input{2,1} = region;
                expres = ['i1.*(i2>0)'];
                Pfname_out = fullfile(dname,[fnames{f,1} '_' peaks_nme{1,p} '_' cond_nme{1,c} '_' region_nme '.nii']);
                if ~exist(Pfname_out,'file') 
                    Output = spm_imcalc_ui(input,Pfname_out,expres);
                end
            end
        end
    end
end

Nfiles = length(fnames);

results = cell(Nfiles,no_peaks,no_cond,Nreg);

for r = 1:Nreg
    region = fullfile(Rdata,regions{r,1});
    Rnii = load_nii(region);
    Rsize = ((sqrt(length(find(Rnii.img>0)))/2)^2); 
    for f = 1:length(fnames)
        for p = 1:no_peaks
            for c = 1:no_cond
                Pfname = fullfile(dname,[fnames{f,1} '_' peaks_nme{1,p} '_' cond_nme{1,c} '_' region_nme]);
                Pdata = [Pfname '.nii'];
                nii = load_nii(Pdata);
                nii_all = nii.img(nii.img~=0);
                nii_mean = sum(nii_all)/Rsize;
                results{f,p,c,r} = nii_mean;
            end
        end
    end
end

results = squeeze(results);
results = permute(results,[2 1 3]);
results2D = reshape(results,f*c,r);
save ROIs.mat results2D;
