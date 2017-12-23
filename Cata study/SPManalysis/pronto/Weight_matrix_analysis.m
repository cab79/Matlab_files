function Weight_matrix_analysis
dbstop if error
close all
%% general inputs for all analyses
% directory for all prt analyses:
D.stats_path = 'C:\Data\Catastrophising study\SPMstats\pronto';
% analysis folder prefix
D.pref = 't-3000_-2_b-3000_-2500_m_-2500_-1000_Grp_Exp_Subject_orig_cleaned_SPNall_prt_';
%D.pref = 't-3000_-2_b-3000_-2500_m_-2500_-1000HighExp_Exp_effect_orig_cleaned_SPNall_prt_gpr_MC';
% analysis folder suffix
%D.suff = {'grpHighExp_gpc_ROI','Exp_gpc_ROI','ExpLL_gpc_ROI'};
D.suff = {'Gender1_Exp_gpc_ROI','Gender2_Exp_gpc_ROI'};
%D.suff = {''};
D.maxarray = 1e+09;

%% get filenames (weights and PRT)
D = get_fnames(D)

%% get image file names
% second input: indices of D.suff to perform this function on
%D = get_image_list(D,1)

%% create projections
% second input: indices of D.suff to perform this function on
%create_projections(D,1)

%% regression of weight vectors
% input 2: indices of D.suff for dependent variable
% input 3: indices of D.suff for independent variable(s)
use_bootstrapped=1;
regress_images(D,{'wimg','pwimg'},1,2,use_bootstrapped)

end

function D = get_fnames(D)

for i = 1:length(D.suff)
    % get weight image path/name
    f = dir(fullfile(D.stats_path,[D.pref D.suff{i}], 'weights*.img'));
    D.wimg{1,i} = fullfile(D.stats_path,[D.pref D.suff{i}], f.name);
    % get weight image path/name
    D.pwimg{1,i} = fullfile(D.stats_path,[D.pref D.suff{i}], ['proj_' f.name]);
    % get prt path/name
    D.prt{1,i} = fullfile(D.stats_path,[D.pref D.suff{i}], 'PRT.mat');
end
end

function D = get_image_list(D,ind)

for i = ind
    % get ERP image path/name per subject/condition from PRT.mat
    load(D.prt{i});
    si=0;
    Ngrp = length(PRT.group);
    for g = 1:Ngrp
        Nsub = length(PRT.group(g).subject);
        for s = 1:Nsub
            Nscans = size(PRT.group(g).subject(s).modality.scans,1);
            for sc = 1:Nscans
                si=si+1;
                D.img{i}{si,1} = PRT.group(g).subject(s).modality.scans(sc,1:end-2);
            end
        end
    end
end
end


function create_projections(D,ind) % UNFINISHED: cannot calculate cov of large arrays
% creates weighted images for group data

for i = ind
    % load weight image and vectorise
    W=spm_vol(D.wimg{i});
    Wimg=spm_read_vols(W(end));
    vWimg = Wimg(:);
    nonan = ~isnan(vWimg);
    vWimg = vWimg(nonan);
        
    % get feature set
    load(D.prt{i});
    fs=PRT.fas.dat(:,:);
    
    % get indices of feature set
    %fsind = PRT.fas.idfeat_img;
    
    % put feature set in weight image dimensions and convert to 2D
    fsd = repmat(double(~isnan(Wimg)),1,1,1,size(fs,1));
    fsd=permute(fsd,[4 1 2 3]);
    fsd(find(fsd)) = fs(:);
    fsd = reshape(fsd,size(fs,1),[]);
    fsd=fsd(:,nonan);
    
    n = size(fsd,1);
    fsd = bsxfun(@minus,fsd,sum(fsd,1)/n);
    % apply weight matrix to img covariance
    if size(fsd,2)>3e4
        % Slower way of estimating but using less memory.
        % size limit setup such that no more than ~1Gb of mem is required:
        % 1Gb/3(nr of matrices)/8(double)= ~40e6 -> sqrt -> 6e3 element
        vwImg = nan(size(vWimg));
        
        chunksize = round(D.maxarray/size(fsd,2));
        nchunk = ceil(size(fsd,2)/chunksize);
        st = 1;
        chunk=[];
        for ic=1:nchunk
            st = st+size(chunk,2);
            en = min(st+chunksize-1,size(fsd,2));
            chunk = fsd(:,st:en);
            %use_col = find(all(~isnan(chunk)));
            %if isempty(use_col)
            %    continue
            %end
            %covmat = nan(size(fsd,2),size(chunk,2));
            %covmat(:,use_col) = (fsd.' * chunk(:,use_col)) ./ (n - 1);
            covmat = (fsd.' * chunk) ./ (n - 1);
            
            %for ir=1:size(fsd,2)
                %covmat = cov(fsd(:,ic),fsd(:,ir),'partialrows');
                %covrow(1,ir) = covmat(1,2);
            %end
            
            vwImg(st:en,1) = covmat'*vWimg; %Phi(:,ic) + kern_vols' * kern_vols(:,ic);
            
            %if any(~isnan(vwImg(st:en,1)))
            %    nn = 'no nans!';
            %else
            %    nn = '';
            %end
            
            disp(['ic = ' num2str(ic*100/nchunk) '%'])
        end
    else
        vwImg = cov(fsd)*vWimg; 
    end
    
    % normalise
    % https://arxiv.org/ftp/arxiv/papers/1606/1606.02840.pdf
    vwImg = vwImg*inv(cov(fsd*vwImg));
    
    pImg = Wimg;
    pImg(nonan)=vwImg;
    
    % write projection image
    Vw = W(1);
    [pth nme ext] = fileparts(D.wimg{i});
    Vw.fname = fullfile(pth,['proj_' nme ext]);
    spm_write_vol(Vw,pImg);
        
end
end


function regress_images(D,fields,indDV,indIV,boot)

R = struct();
for f = 1:length(fields)
    R.DV = D.(fields{f}){indDV};
    if ~exist(R.DV)
        disp('no projections');
        continue
    end
    % load DV
    V=spm_vol(R.DV);
    Vdv=spm_read_vols(V(end));
    % vectorise
    dv = Vdv(:);
    
    for i = 1:length(indIV)
        % load IV
        R.IV{i} = D.(fields{f}){indIV(i)};
        V=spm_vol(R.IV{i});
        Viv=spm_read_vols(V(end));
        % vectorise
        iv(:,i) = Viv(:);
    end
    dvi = ~isnan(dv);
    ivi = ~isnan(iv(:,1));
    if sum(dvi)~=sum(ivi)
        error('dv and iv are different sizes')
    end
    % multiple regression
    %STATS contains:
    %the R-square statistic, 
    %the F statistic 
    %the p value
    %an estimate of the error variance.
    [b,bint,r,rint,R.st.stats] = regress(dv(dvi),[ones(sum(ivi),1) iv(ivi,:)]);
    if boot && strcmp(fields{f},'wimg')
        DVs = dir(fullfile(D.stats_path,[D.pref D.suff{indDV}],'perm_weights*'));
        DV = fullfile(D.stats_path,[D.pref D.suff{indDV}],DVs.name);
        for ii = 1:length(indIV)
            IVs = dir(fullfile(D.stats_path,[D.pref D.suff{indIV(ii)}],'perm_weights*'));
            IV{ii,1} = fullfile(D.stats_path,[D.pref D.suff{indIV(ii)}],IVs.name);
        end
        P = regress_permuted(DV,IV);
        R.st.permresults = extractfield(P, 'rspacetime'); 
        R.st.pval = sum(abs(R.st.permresults) > abs(R.st.stats(1))) / length(R.st.permresults);
        pID='permuted';
        R.st.CIs = prctile(R.st.permresults,[2.5 97.5]);
    else
        R.st.pval = R.st.stats(3);
        pID='';
        R.st.CIs = [nan,nan];
    end
    figure;scatter(iv(ivi,:),dv(dvi))
    title([fields{f} ': space-time r2=' num2str(R.st.stats(1)) ', ' pID ' p=' num2str(R.st.pval) ', CIs = ' num2str(R.st.CIs(1)) ' ' num2str(R.st.CIs(2))])
    
    % regression over space for each time point
    R.s.rtime=nan(size(Vdv,3),1);
    R.s.pval = nan(size(Vdv,3),1);
    R.s.CIs = nan(size(Vdv,3),2);
    for d = 1:size(Vdv,3)
        if all(isnan(Vdv(:,:,d)))
            continue;
        end
        dvt = Vdv(:,:,d);
        dvt = dvt(~isnan(dvt));
        for i = 1:size(iv,2)
            ivtemp = Viv(:,:,d);
            ivt(:,i) = ivtemp(~isnan(ivtemp));
        end
        [b,bint,r,rint,stats] = regress(dvt,[ones(length(dvt),1) ivt]);
        R.s.rtime(d) = stats(1);
        if boot && strcmp(fields{f},'wimg')
            permresults = [];
            for i = 1:length(P)
                permresults(i) = P(i).rtime(d);
            end
            R.s.pval(d) = sum(abs(permresults) > abs(R.s.rtime(d))) / length(permresults);
            pT='permuted';
            R.s.CIs(d,:) = prctile(permresults,[2.5 97.5]);
        else
            R.s.pval = stats(3);
            pT='';
            R.s.CIs(d,:) = [nan,nan];
        end
    end
    figure;
    rtime = R.s.rtime;
    rtime_corr = nan(size(rtime));
    if strcmp(pT,'permuted')
        rtime(R.s.pval>0.05)=nan;
        [R.s.pID,R.s.pN] = fdr(R.s.pval,0.05);
        if ~isempty(R.s.pN)
            rtime_corr = rtime;
            rtime_corr(rtime_corr>R.s.pN) = nan;
        end
    end
    plot(rtime,'b'); hold on
    plot(rtime_corr,'r'); hold off
    title([fields{f} ': r2 over time, ' pID])
    [maxr,maxt]=max(rtime);
    xx=Viv(:,:,maxt);
    yy=Vdv(:,:,maxt);
    figure;scatter(xx(:),yy(:));
    title([fields{f} ': max r2 @dp ' num2str(maxt) ' =' num2str(maxr)])
    
    % weights per time ROI
    % get weights from DV and IV prt
    load(D.prt{indDV});
    dvwt = PRT.model.output.weight_ROI{1, 1}(:,end);
    for i = 1:length(indIV)
        load(D.prt{indIV});
        ivwt(:,i) = PRT.model.output.weight_ROI{1, 1}(:,end);
        if length(ivwt(:,i))~=length(dvwt)
            error('dv and iv are different sizes')
        end
    end
    % multiple regression
    %STATS contains:
    %the R-square statistic, 
    %the F statistic 
    %the p value
    %an estimate of the error variance.
    [b,bint,r,rint,R.t.stats] = regress(dvwt,[ones(length(dvwt),1) ivwt]);
    figure;scatter(ivwt,dvwt)
    title([fields{f} ': weights over time, r2=' num2str(R.t.stats(1)) ', p=' num2str(R.t.stats(3))])
end
save(fullfile(D.stats_path,[D.pref D.suff{indDV}],'weights_stats.mat'),'R');
end

function P = regress_permuted(DV,IV)
    
% DV & IV: path of permuted weights folder
dimgs = dir(fullfile(DV,'*.img'));
for i = 1:length(dimgs)
    disp(['regress_permute: loading DV weight image ' num2str(i)])
    DVimg = fullfile(DV,dimgs(i).name);
    V=spm_vol(DVimg);
    P(i).Vdv=spm_read_vols(V);
    P(i).dv = P(i).Vdv(:);
    P(i).dvi = ~isnan(P(i).dv);
end
for ii = 1:length(IV)
    iimgs = dir(fullfile(IV{ii},'*.img'));
    if length(dimgs)~=length(iimgs)
        error('different number of permutations for dv and iv')
    end
    for i = 1:length(iimgs)
        disp(['regress_permute: loading IV' num2str(ii) ' weight image ' num2str(i)])
        IVimg = fullfile(IV{ii},iimgs(i).name);
        V=spm_vol(IVimg);
        P(i).Viv{ii}=spm_read_vols(V);
        P(i).iv(:,ii) = P(i).Viv{ii}(:);
        P(i).ivi(:,ii) = ~isnan(P(i).iv(:,ii));
        if sum(P(i).dvi)~=sum(P(i).ivi(:,ii))
            error('dv and iv are different sizes')
        end
    end
end
for i = 1:length(dimgs)
    % multiple regression
    %STATS contains:
    %the R-square statistic, 
    %the F statistic 
    %the p value
    %an estimate of the error variance.
    [b,bint,r,rint,stats] = regress(P(i).dv(P(i).dvi),[ones(sum(P(i).ivi(:,1)),1) P(i).iv(P(i).ivi(:,1))]);
    P(i).rspacetime = stats(1);
end

% regression over space for each time point
for i = 1:length(dimgs)
    disp(['regress_permute: regression over space for weight image ' num2str(i)])
    P(i).rtime=nan(1,size(P(i).Vdv,3));
    for d = 1:size(P(i).Vdv,3)
        if all(isnan(P(i).Vdv(:,:,d)))
            continue;
        end
        dvt = P(i).Vdv(:,:,d);
        dvt = dvt(~isnan(dvt));
        for ii = 1:size(P(i).iv,2)
            ivtemp = P(i).Viv{ii}(:,:,d);
            ivt(:,ii) = ivtemp(~isnan(ivtemp));
        end
        [b,bint,r,rint,stats] = regress(dvt,[ones(length(dvt),1) ivt]);
        P(i).rtime(d) = stats(1);
    end
end
end