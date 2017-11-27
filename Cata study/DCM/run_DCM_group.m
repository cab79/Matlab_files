function run_DCM_group(S,DCM)

%% RUN
cd(S.outpath)
%setup save name and run info
sii=1;
if exist(['Run_info' num2str(S.run_num) '.mat'],'file')
    subs=[];
    load(['Run_info' num2str(S.run_num) '.mat']);
    for si = 1:length(S.run_subject_subsets)
        for su = 1:length(subs)
            ei(su) = isequal(subs(su),S.run_subject_subsets(si));
        end
        if ~any(ei)
            sii=si;
            break
        end
    end
    subs(1,length(subs)+1)=S.run_subject_subsets(sii);
else 
    subs=S.run_subject_subsets(sii);
end
save(['Run_info' num2str(S.run_num) '.mat'],'subs');
subjects = S.run_subject_subsets{sii}

if length(S.run_subject_subsets)>1
    S.sub_ind = S.run_subject_subsets{sii};
end

% options
if size(S.design,2)>1 && size(S.design,1)==1
    S.design = S.design';
end
DCM.options.Tdcm(1)  = S.timewin(1);     % start of peri-stimulus time to be modelled
DCM.options.Tdcm(2)  = S.timewin(2);   % end of peri-stimulus time to be modelled
DCM.options.trials   = unique(S.conds(S.conds>0));
DCM.xU.X = S.design;
DCM.xU.name = {S.contrastname};

%--------------------------------------------------------------------------
% Specify connectivity model
%--------------------------------------------------------------------------

% A matrix is the connections of the first trial/event
% NB: column index corresponds to the source area, and the row index to the
% target area.

% forward connections
DCM.A{1} = xlsread(S.model_file,'A_forward');
% backward connections
DCM.A{2} = xlsread(S.model_file,'A_backward');
% lateral connections
DCM.A{3} = xlsread(S.model_file,'A_lateral');
% B matrix: Gain modulations of connection strengths as set in the A-matrices
% Models the difference between the first and the other modelled evoked responses.
% For example, for two evoked responses, DCM explains the first response by
% using the A-matrix only. The 2nd response is modelled by modulating these connections 
% by the weights in the B-matrix.
% Set a single B model here to only test one model. 
%DCM.B{1} = DCM.A{1} + DCM.A{2};
%DCM.B{1}(1,1) = 1;
%DCM.B{1}(2,2) = 1;
DCM.B{1} = xlsread(S.model_file,'B');
% Alternatively, to test any combination of forward, backward, lateral and intrinsic connections in
% multipe models, type DCM.Bc = [1,2,3,4] for all, or for a subset. Leave
% empty to use the singel model defined above.
DCM.Ac = xlsread(S.model_file,'A_Acomb'); % order: F,B,L (no I)
DCM.Ac = find(DCM.Ac);
DCM.Anc = 1:3;
DCM.Anc(DCM.Ac) = [];
DCM.Bc = xlsread(S.model_file,'B_Acomb'); % order: F,B,L,I
DCM.Bc = find(DCM.Bc);
DCM.Bnc = 1:4;
DCM.Bnc(DCM.Bc) = [];
% C matrix: Inputs - can be to one or many areas
Cval = xlsread(S.model_file,'C');
DCM.C = Cval(1,:)';
% indices of nodes that should be tested for self-connection changes
DCM.Bi = Cval(2,:); 
DCM.Bi = find(DCM.Bi);

% read in locations data
[Lpos,Sname]=xlsread(fullfile(S.outpath,S.loc_file));
% specify the prior source locations (in mm in MNI coordinates). load the prior locations from a file
DCM.Lpos  = Lpos';
% enter the source names (one name in one row). 
DCM.Sname = Sname';
Nareas    = size(DCM.Lpos,2);

% NOTATION FOR SOME KEY FIELD NAMES:
%     DCM{i}.M.pE	- prior expectation of parameters
%     DCM{i}.M.pC	- prior covariances of parameters
%     DCM{i}.Ep		- posterior expectations
%     DCM{i}.Cp		- posterior covariance
%     DCM{i}.F		- free energy

if ~exist(S.outpath,'dir')
    mkdir(S.outpath);
end
cd(S.outpath)

if length(S.run_subject_subsets)>1
    sname = ['GCM_fit' num2str(S.run_num) num2str(S.run_subject_subsets{sii}(1)) '_' num2str(S.run_subject_subsets{sii}(end))];
else
    sname = ['GCM_fit' num2str(S.run_num)];
end

if S.loadorsave && exist(fullfile(S.outpath,[sname '.mat']),'file')
    %load(fullfile(S.outpath,[sname '.mat']));
else
    [~,~,pdata] = xlsread(S.pdatfile);
    grp_col = find(strcmp(pdata(1,:),S.grphead));
    sub_col = find(strcmp(pdata(1,:),S.subhead));
    inc_col = find(strcmp(pdata(1,:),S.inchead));

    %Change char Grp inputs to numbers
    grpdat = pdata(2:end,grp_col);
    if isnumeric(grpdat{2,1})
        grptype = unique([grpdat{:}]);
        grptype(grptype==0)=[];
        Ngrp = length(grptype);
    else
        grptype = unique(grpdat);
        grptype(isempty(grptype))=[];
        Ngrp = length(grptype);
        for g = 1:Ngrp
            grp_idx = cellfun(@(x) any(strcmp(grptype(g),x)), grpdat, 'UniformOutput', 0);
            grpdat(cell2mat(grp_idx)) = {[g]};
        end
    end

    % find index of subjects to include in analysis
    SubInd = cell(Ngrp,1);
    Subs = [];
    inc_idx = cellfun(@(x) ismember(x,S.include_codes), pdata(2:end,inc_col), 'UniformOutput', 0);
    inc_idx = find(cell2mat(inc_idx));

    % find subject indices for each specific group
    for g = 1:Ngrp
        grp_idx = find(cellfun(@(x) x==g, grpdat, 'UniformOutput', 1));
        SubInd{g,1} = intersect(inc_idx,grp_idx);
        Nsub(g,1) = length(SubInd{g,1});
    end

    % create SubInd vector
    SubInd_ana=cell2mat(SubInd);
    % reduce subjects to those selected in this batch
    if S.sub_ind
        SubInd_ana = SubInd_ana(S.sub_ind);
    end

    % between subject effects: constant, group difference
    %--------------------------------------------------------------------------
    if size(S.grps,1)==1
        Xb = [ones(sum(Nsub),1) ones(sum(Nsub),1)];
    elseif size(S.grps,1)==2
        Xb = [ones(sum(Nsub),1) [ones(Nsub(1),1);-ones(Nsub(2),1)]];
    end

    % covariates
    if ~isempty(S.cov_names)
        if ~strcmp(S.cov_names{1},'')
            for c = 1:length(S.cov_names)
                cov_col = find(strcmp(pdata(1,:),S.cov_names{c}));
                covdat = pdata(2:end,cov_col);

                % convert categorical (char) covariate data to numbers
                if ~isnumeric(covdat{2,1})
                    covtype = unique(covdat);
                    Ncovtype = length(covtype);
                    for c = 1:Ncovtype
                        cov_idx = cellfun(@(x) any(strcmp(covtype(c),x)), covdat, 'UniformOutput', 0);
                        covdat(cell2mat(cov_idx)) = {[c]};
                    end
                end
                covdat = covdat(vertcat(SubInd_ana{:}));

                Xb = horzcat(Xb,[covdat{:}]');
            end
        end
    end
    
    
    Cv = unique(DCM.C); % c values
    Ci = find(double(~(Cv==0)) .* double(Cv<100)); % c indices
    nC = length(Ci); % number of values
    C=[];
    Cr=[];
    if nC>1
        k = spm_perm_mtx(nC);
        % for each possible combination of factors in A
        for i = 1:(2^nC);
            C{i}     = sparse(1,Nareas); % sparse 5x5 double array - for connecting 5 nodes
            Cr{i} = C{i};
            % for each factor in C
            for f = 1:nC
                try
                    if k(i,f)
                        C{i} = C{i} + double(DCM.C==Cv(Ci(f)))' + double(DCM.C==100)'; 
                        Cr{i} = Cr{i} + (DCM.C.*double(DCM.C==Cv(Ci(f))))'; 
                    end
                end
            end
            % ensure all priors equal 1
            C{i} = C{i}>0;
            % remove model if no inputs
            if sum(C{i})==0
                C(i)=[];
            end
        end
        % remove models with too many connections
        delC=zeros(1,length(C));
        sizCfm = length(find(~(DCM.C==0)));
        for i = 1:length(C);
            % identify the full model to keep
            sizC=length(find(~(C{i}==0)));
            if sizC==sizCfm
                Cfm = C{i};
            end
            for f = 1:nC
                if any(Cr{i}==Cv(Ci(f))) && any(Cr{i}==-Cv(Ci(f)))
                    delC(i)=1;
                end
            end
        end
        C(find(delC))=[];
        nC = length(C);
    else
        C{1} = DCM.C';
        nC=1;
        Cfm = C{1};
    end



    B=[];
    if any(DCM.Bc)
        % model space - within subject effects - creates B matrix of models to
        % compare
        % B matrix: Gain modulations of connection strengths as set in the A-matrices
        % In this case, we are creating 8 models: see Table 1 of http://www.sciencedirect.com/science/article/pii/S105381191501037X#s0015
        %----------------------------------------------------------------------
        % create (2^n x n) sparse matrix of indices permuted over n
        k = spm_perm_mtx(length(DCM.Bc));
        % for each reduced model, create sparse matrices of ones
        % sparse matrices are used when we expect most elements to contain zeros, and only a few non-zero elements.
        % this reduces memory size. In this case, it's done as an efficient way of
        % coding the connections (values of ones)
        for i = 1:(2^length(DCM.Bc));
            B{i}     = sparse(Nareas,Nareas); % sparse 5x5 double array - for connecting 5 nodes
            % B comprises every combination (8 in total) of the options below: 
            try
                if k(i,1) && any(DCM.Bc==1)
                    B{i} = B{i} + DCM.A{1}; % All forward connections
                end
            end
            try
                if k(i,2) && any(DCM.Bc==2)
                    B{i} = B{i} + DCM.A{2}; % All backward connections
                end
            end
            try
                if k(i,3) && any(DCM.Bc==3)
                    B{i} = B{i} + DCM.A{3}; % All lateral connections
                end
            end
            try
                if k(i,4) && any(DCM.Bc==4)
                    B{i} = B{i} + sparse(DCM.Bi,DCM.Bi,1,Nareas,Nareas); % All intrinsic connections
                end
            end
            % ensure all priors equal 1
            B{i} = B{i}>0;
            % add the zero elements in
            B{i}     = full(B{i});
        end
        Bfm=B(1);
        nB = length(B);
    elseif any(DCM.B{1})
        Bv = unique([DCM.B{:}]); % B values
        Bi = find(double(~(Bv==0)) .* double(Bv<100)); % B indices
        nB = length(Bi); % number of values
        B=[];
        Br=[];
        k = spm_perm_mtx(nB);
        % for each possible combination of factors in B
        for i = 1:(2^nB);
            B{i}     = sparse(Nareas,Nareas); % sparse 5x5 double array - for connecting 5 nodes
            Br{i}     = B{i}; 
            % for each factor in B
            for f = 1:nB
                try
                    if k(i,f)
                        B{i} = B{i} + double([DCM.B{:}]==Bv(Bi(f))) + double([DCM.B{:}]==-1); 
                        Br{i} = Br{i} + [DCM.B{:}].*double([DCM.B{:}]==Bv(Bi(f))); 
                    end
                end
            end
            % ensure all priors equal 1
            B{i} = B{i}>0;
            % remove any with no connections
            if sum(B{i})==0
                B(i)=[];
            end
        end
        % remove models with too many connections
        delB=zeros(1,length(B));
        sizBfm = length(find(~([DCM.B{:}]==0)));
        for i = 1:length(B);
            % identify the full model to keep
            sizB=length(find(~([B{i}]==0)));
            if sizB==sizBfm
                Bfm = B{i};
            end
            for f = 1:nB
                if any(any(Br{i}==Bv(Bi(f)))) && any(any(Br{i}==-Bv(Bi(f))))
                    delB(i)=1;
                end
            end
        end
        B(find(delB))=[];
        
        % remove models in which connections exist in higher levels of the
        % hierarchy but are not connected to regions in which there are
        % inputs.
        delB=zeros(1,length(B));
        inp = find(Cfm);
        for i = 1:length(B);
            % for each column in the connection matrix except input
            % cols
            cols = 1:size(B{i});
            cols(inp) = [];
            for col = cols
                % if there are any forward connections
                if any(B{i}(:,col))
                    % while those regions are not yet found connected to
                    % the input region
                    allconrow=[];
                    newcol=inp; % initialise newcol to input regions
                    while ~any(ismember(allconrow,col))
                        conrow = [];
                        for cc = 1:length(newcol)
                            % find all connected rows
                            conrow = [conrow;find(B{i}(:,newcol(cc)))];
                        end
                        conrow = unique(conrow);
                        % if there are no connections found, break
                        if isempty(conrow) || all(ismember(conrow,allconrow))
                            break
                        else
                            allconrow = [allconrow;conrow];
                            newcol = conrow;
                        end
                    end
                    % mark model for deletion if this col is not
                    % connected to the input region
                    if ~any(ismember(allconrow,col))
                        delB(i)=1;
                        break
                    end
                end
            end
        end
        B(find(delB))=[];
        nB = length(B);
    else
        nB=1;
        B{1} = sparse(Nareas,Nareas);
        Bfm = B;
    end

    %model space for A - must be no B matrix to use this
    %nA = sum(unique([DCM.A{:}])>0);
    Av = unique([DCM.A{DCM.Ac}]); % A values
    Ai = find(double(~(Av==0)) .* double(Av<100)); % A indices
    nA = length(Ai); % number of values
    A=[];
    Ar=[];
    if sum(DCM.B{:}(:))==0 && nA>1
        k = spm_perm_mtx(nA);
        % for each possible combination of factors in A
        for i = 1:(2^nA);
            %for each connection type in A
            for c = DCM.Ac
                A{i}{c} = sparse(Nareas,Nareas); % sparse 5x5 double array - for connecting 5 nodes
                Ar{i}{c} = A{i}{c}; % real valued version of A
                % for each factor in A
                for f = 1:nA
                    try
                        if k(i,f)
                            A{i}{c} = A{i}{c} + double(DCM.A{c}==Av(Ai(f))) + double(DCM.A{c}==100); 
                            Ar{i}{c} = Ar{i}{c} + DCM.A{c}.*double(DCM.A{c}==Av(Ai(f))); 
                        end
                    end
                end
                % ensure all priors equal 1
                A{i}{c} = A{i}{c}>0;
            end
            for c = DCM.Anc
                A{i}{c} = DCM.A{c}; 
                % ensure all priors equal 1
                A{i}{c} = A{i}{c}>0;
            end
            % remove any with no connections
            if sum(A{i}{c})==0
                A(i)=[];
            end
        end
        
        % remove all other models with too many connections
        % i.e. both positive and negative values of the same number, e.g. 2
        % AND -2 - these models should be eliminated
        delA=zeros(1,length(A));
        sizAfm = length(find(~([DCM.A{DCM.Ac}]==0)));
        for i = 1:length(A);
            % identify the full model to keep
            sizA=length(find(~([A{i}{DCM.Ac}]==0)));
            if sizA==sizAfm
                Afm = A{i};
            end
            %for each connection type in A
            for c = DCM.Ac
                for f = 1:nA
                    if any(any(Ar{i}{c}==Av(Ai(f)))) && any(any(Ar{i}{c}==-Av(Ai(f))))
                        delA(i)=1;
                    end
                end
            end
        end
        A(find(delA))=[];
        
        % remove models in which connections exist in higher levels of the
        % hierarchy but are not connected to regions in which there are
        % inputs.
        delA=zeros(1,length(A));
        inp = find(Cfm);
        for i = 1:length(A);
            disp(['checking connections of model ' num2str(i)]);
            %for forward connections only
            for c = 1
                % for each column in the connection matrix except input
                % cols
                cols = 1:size(A{i}{c});
                cols(inp) = [];
                for col = cols
                    % if there are any connections
                    if any(A{i}{c}(:,col))
                        disp(['column ' num2str(col)]);
                        % while those regions are not yet found connected to
                        % the input region
                        allconrow=[];
                        newcol=inp; % initialise newcol to input regions
                        while ~any(ismember(allconrow,col))
                            conrow = [];
                            for cc = 1:length(newcol)
                                % find all connected rows
                                conrow = [conrow;find(A{i}{c}(:,newcol(cc)))];
                            end
                            conrow = unique(conrow);
                            % if there are no new connections found, break
                            if isempty(conrow) || all(ismember(conrow,allconrow))
                                break
                            else
                                allconrow = [allconrow;conrow];
                                newcol = conrow;
                            end
                        end
                        % mark model for deletion if this col is not
                        % connected to the input region
                        if ~any(ismember(allconrow,col))
                            delA(i)=1;
                            break
                        end
                    end
                end
            end
        end
        A(find(delA))=[];
        nA = length(A);
    else
        A{1}{1} = sparse(Nareas,Nareas);
        A{1}{2} = sparse(Nareas,Nareas);
        A{1}{3} = sparse(Nareas,Nareas);
        nA=1;
        Afm = A{1};
    end


    % combinations of reduced models of A, B and C
    [Aind,Bind,Cind] = meshgrid(1:nA, 1:nB, 1:nC);
    combs = [Aind(:),Bind(:),Cind(:)];
    nMod = length(combs);

    % loop through subjects to set up for inversion
    %--------------------------------------------------------------------------
    par=[];
    if S.prepare_data==0 && exist('DCM_PEB_alldata.mat','file')
        load DCM_PEB_alldata
    else
        subID={};
        files={};
        GCM=[];
        i=0;
        load(fullfile(S.fid_dir,'meegfid2.mat'));
        for i = 1:length(SubInd_ana)
            subID = pdata{SubInd_ana(i)+1,sub_col};
            if isnumeric(subID); subID = num2str(subID); end;
            fname = dir(fullfile(S.filepath,[S.fpref '*' subID '*' S.fmid  '*' S.fsuff]));
            files = fname.name;
            % Baseline Correction
            clear SB
            SB.D = fullfile(S.filepath,files);
            SB.timewin = S.basewin;
            SB.save = 0; % save in separate file
            SB.prefix = 'b'; % for output, only if saving a new file
            spm_eeg_bc(SB);

            % Data filename
            if isfield(DCM,'xY'); DCM = rmfield(DCM,'xY');end
            DCM.xY.Dfile = fullfile(S.filepath,files);

            % prepare file: order conditions
            %if isfield(S,'sortconds')
            %    if ~isempty(S.sortconds)
            %        matlabbatch{1}.spm.meeg.preproc.prepare.D = {DCM.xY.Dfile};
            %        matlabbatch{1}.spm.meeg.preproc.prepare.task{1}.sortconditions.label = S.sortconds;
            %        spm_jobman('initcfg')
            %        spm_jobman('run',matlabbatch);
            %    end
            %end

            % add fiducials
            load(DCM.xY.Dfile)
            if ~isequal(D.fiducials,meegfid)
                D.fiducials = meegfid;
                save(DCM.xY.Dfile,'D');
            end

            for j = 1:nMod+1 % all combinations of reduced models, plus full model
                disp(['preparing data - subject number: ' subID ', model:' num2str(j)]); % display progress 
                if isfield(DCM,'M'); DCM = rmfield(DCM,'M');end
                DCM.name = ['DCM_' DCM.options.analysis '-ana_' DCM.options.model '-model_' DCM.options.spatial '-spatial_' num2str(j) '-model'];

                if j==1
                    % Data and spatial model
                    DCM  = spm_dcm_erp_data_CAB(DCM,1,S);
                    xY =DCM.xY;
                    % full model
                    DCM.A   = Afm;
                    DCM.B   = Bfm;
                    DCM.C   = Cfm';
                else
                    DCM.xY =xY;
                    % reduced models
                    DCM.A   = A{combs(j-1,1)};
                    DCM.B   = B(combs(j-1,2));
                    DCM.C   = C{combs(j-1,3)}';
                end
                GCM{i,j} = DCM;

                % find number of free parameters to determine which is the full model
                %[ii,rC,rE,Np] = spm_find_pC(GCM{i,j});
                %par(:,j)     = sparse(ii,1,true,Np,1);
            end
            %[~,fu] = max(sum(par));
            fu=1;
            if i==1
                % Spatial model
                GCM{i,fu} = spm_dcm_erp_dipfit(GCM{i,fu},1);
                M.dipfit = GCM{i,fu}.M.dipfit;
            else
                GCM{i,fu}.M = M;
            end
        end
        if S.save_data
            disp('Saving data - can take a while...'); % display progress    
            save DCM_PEB_alldata GCM -v7.3
        end
    end

    % invert rreduced models (standard inversion) - takes HOURS!!!
    %==========================================================================


    %----------------------------------------------------------------------


    if S.invert_all
        if any(S.run_model_subsets>0)
            GCM_temp = GCM(:,S.run_model_subsets);
            GCM_temp = spm_dcm_fit(GCM_temp); % calls spm_dcm_erp for each subject and model
            GCM(:,S.run_model_subsets) = GCM_temp;
        else
            GCM = spm_dcm_fit(GCM); % calls spm_dcm_erp for each subject and model
        end
    else
        %if length(S.run_subject_subsets)>1
        %    GCM_temp = spm_dcm_fit(GCM(S.run_subject_subsets{sii},fu)); % calls spm_dcm_erp for each subject and model
        %    GCM(S.run_subject_subsets{sii},fu) = GCM_temp;
        %else
            GCM(:,fu) = spm_dcm_fit(GCM(:,fu),S); % calls spm_dcm_erp for each subject and model
        %end
    end
    subjects = S.run_subject_subsets{sii};
    excl_full=S.excl_full;
    save([sname '.mat'],'GCM','DCM','Xb','subjects','excl_full','-v7.3');
    save(['models_' num2str(S.run_num) '.mat'],'A','B','C','combs','-v7.3');
end

if S.run_PEB
    lname = ['GCM_fit' num2str(S.run_num)];
    files = dir([lname '*.mat']);
    if length(files)>1
        if exist([lname '.mat'],'file')
            load([lname '.mat']);
            GCM_all=GCM;
            st=2;
        else
            st=1;
        end
        for f = st:length(files)
            load(files(f).name)
            GCM_all(subjects,:) = GCM;
            clear GCM
        end
        GCM=GCM_all;
        clear GCM_all
        save([sname '.mat'],'GCM','DCM','Xb','subjects','A','B','C','combs','-v7.3');
        for f = st:length(files)
            delete(files(f).name)
        end
    end
    fu=1;
    clear GCM DCM % to free memory
    DCM_PEB_trim([lname '.mat'],S.outpath,S.run_num,fu,S.loadorsave)
end
