clear all
%close all
Pbase = 'F:\Dell\bloodA\Image_analysis_files\examples';
PsubPET = 'PET';
scans = {1 2 'diff'};
scan_name = {'pain','nopain'};

Ebase = 'C:\Data\PET-LEP\Preprocessed';
[NUM,TXT,RAW] = xlsread(fullfile(Ebase,'gfp_results_ERP_evoked_1_2_eachcond_nonorm.xlsx'));
peak_no = 24;
readcol=peak_no+2;

run('C:\Data\Matlab\Matlab_files\PET-LEP\EEG scripts\loadsubj.m');
grplist = [2 3,... % healthy
            5 6];   % patient
select_scan= [1]; %(1=pain, 2=nopain)
subjects = subjlists(grplist);
        
scanlist=cell(1,1);
Nsub=0;


%% PRONTO
PRT = Design_batch_pronto(D)

load(fullfile(D.batch,'PRT.mat'));

for s = 1:length(subjects)
     for s2 = 1:size(subjects{s,1},1) 
         
        % get EEG data
        if exist('RAW','var')
            EEGnme = strsplit(subjects{s,1}{s2,1},'_');
            EEGnme = EEGnme{1};
            EEGsubidx = find(strcmp(RAW(:,1),EEGnme));
            rt_subj = RAW{EEGsubidx,readcol};
        else
            rt_subj=[];
        end

        %identify and list PET file names
        Pfiles = dir(fullfile(Pbase, subjects{s,1}{s2,2}, PsubPET,['s8wrh*nnls_Vd_mag-310_rsl.img']));
        Pfiles_sub = dir(fullfile(Pbase, subjects{s,1}{s2,2}, PsubPET,['s8wrh*nnls_Vd_mag-310_rsl_SUBTRACTION.img']));
        Pfiles = [Pfiles;Pfiles_sub];

        if ismember(grplist(s),[2 5])
            scan_ana = subjects{s,1}{s2,3};
            if scan_ana ~= select_scan
                continue
            else 
                scan_ana=1;
            end
        elseif ismember(grplist(s),[3 6])
            scan1 = subjects{s,1}{s2,3};
            scan2 = [scans{1:2}];
            scan2(scan1) = [];
            scan_ana = select_scan;
            scan_order = [scan1,scan2];
            Pfiles(1:2) = Pfiles(scan_order);
        end
        
        Nsub=Nsub+1;
        PRT.group.subject(Nsub).subj_name = subjects{s,1}{s2,2};
        PRT.group.subject(Nsub).modality = PRT.group.subject(1).modality;
        if length(select_scan)>1 % i.e. if there is no within-subject contrast
            infiles=cell(1,1);
            for sa=1:length(scan_ana)
                %PRT.group.subject(Nsub).modality(sa).mod_name = scan_name{sa};
                infiles{sa,1}=[fullfile(Pbase, subjects{s,1}{s2,2}, PsubPET, Pfiles(sa).name) ',1'];
                PRT.group.subject(Nsub).modality.design.conds(sa).scans = sa;
            end
            PRT.group.subject(Nsub).modality.scans = char(infiles);
            PRT.group.subject(Nsub).modality.TR = 1;
            PRT.group.subject(Nsub).modality.design.covar = [];
        else
            PRT.group.subject(Nsub).modality.scans = [fullfile(Pbase, subjects{s,1}{s2,2}, PsubPET, Pfiles(scan_ana).name) ',1'];
            PRT.group.subject(Nsub).modality.rt_subj = rt_subj;
            PRT.group.subject(Nsub).modality.covar = [];
        end
        %scanlist{Nsub,1} = [fullfile(Pbase, subjects{s,1}{s2,2}, PsubPET, Pfiles(select_scan).name) ',1'];
     end
end

save PRT PRT

%% SPM
%load job_allcov;
%matlabbatch{1,1}.spm.stats.factorial_design.des.mreg.scans = pain_scans;
%matlabbatch{1,1}.spm.stats.factorial_design.des.fblock.fsuball.specall.scans = both_scans_meannorm;
%matlabbatch{1,1}.spm.stats.factorial_design.des.fblock.fsuball.specall.scans = both_scans;
%matlabbatch{1,1}.spm.stats.factorial_design.des.fblock.fsuball.specall.scans = nonpain_scans;
%matlabbatch{1,1}.spm.stats.factorial_design.des.fblock.fsuball.specall.scans = [nonpain_scans; pain_scans];
%matlabbatch{1,1}.spm.stats.factorial_design.des.mreg.scans = both_scans_meannorm;
%matlabbatch{1,1}.spm.stats.factorial_design.des.mreg.scans = subtracted_scans;
%matlabbatch{1,1}.spm.stats.factorial_design.des.mreg.scans = mean_scans;
%matlabbatch{1,1}.spm.stats.factorial_design.des.mreg.scans = pain_scans;
%matlabbatch{1,1}.spm.stats.factorial_design.des.mreg.scans = nonpain_scans;
%save job_allcov matlabbatch;
