function Normality_test_residuals(S)
dbstop if error
S.filenames = {'Res*.xlsx'};


if isempty(S.spm_dir)
    S.spm_paths = dir(fullfile(S.spmstats_path,'*spm'));
    S.spm_paths = {S.spm_paths(:).name};
else
    S.spm_paths = {S.spm_dir};
end

for sp = 1:length(S.spm_paths)
    S.spm_path = fullfile(S.spmstats_path,S.spm_paths{sp});

    if ~isfield(S,'Fm');
        load(fullfile(S.spm_path,'SPM.mat'));
        S.Fm = SPM.xX.I; % factor matrix
    end

    % load SPM design and list factors
    load(fullfile(S.spm_path,S.batch));
    S.fact = {matlabbatch{1,1}.spm.stats.factorial_design.des.fblock.fac(:).name};

    S.clus_path={};
    alldir=dir(fullfile(S.spm_path,'*_clusters'));
    for c = 1:length(alldir)
        S.clus_path{c} = fullfile(S.spm_path,alldir(c).name);
    end

    for cldir = 1:length(S.clus_path)

        % factor abbreviations
        factors = strsplit(alldir(cldir).name,'_');
        factors = factors(1:end-1);

        % find columns of S.fact for these factors, so we can find the correct
        % column of the factor matrix
        fact_row = find(ismember(S.fact,factors));

        factind = S.Fm(:,fact_row+1);
        [~,~,rowind] = unique(factind,'rows','stable');

        rfiles = dir(fullfile(S.clus_path{cldir},S.filenames{:}));

        h = NaN(length(unique(rowind)),length(rfiles));
        pval=h;
        for r = 1:length(rfiles)

            % load Res excel file 
            rname = rfiles(r).name;
            [datnum,dattxt,datraw] = xlsread(fullfile(S.clus_path{cldir},rname));

            factdat=datnum(:,end);

            for c=1:length(unique(rowind))
                %h(c,fi) = kstest(dat.Data(cells==c));
                [h(c,r),pval(c,r),~] = swtest(factdat(rowind==c), 0.05);
            end

        end

        norm_head = horzcat({rfiles(:).name});
        out = vertcat(norm_head,num2cell(h),num2cell(pval));
        xlswrite(fullfile(S.clus_path{cldir},['normality_tests.xls']),out);

    end
end
       