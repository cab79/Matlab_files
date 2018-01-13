function Weight_matrix_analysis
dbstop if error
close all
%% general inputs for all analyses
% directory for all prt analyses:
D.stats_path = 'C:\Data\Catastrophising study\SPMstats\pronto';
% analysis folder prefix
%D.pref = 'Group\Per condition\t-3000_-2_b-3000_-2500_m_-2500_-1000_Grp_Exp_Subject_orig_cleaned_SPNall_prt_';
%D.pref = 't-3000_-2_b-3000_-2500_m_-2500_-1000HighExp_Exp_effect_orig_cleaned_SPNall_prt_gpr_MC';
%D.pref = 'Group Expectancy\t-5500_1500_b-5500_-5000_m_0_1500_Grp_Exp_Subject_orig_cleaned_trialNmatch_prt';
D.pref = '';
% analysis folder suffix
D.suff = {'t-5500_1500_b-5500_-5000_m_0_1500_Grp_Exp_Subject_orig_cleaned_trialNmatch_prt_Int_gpc_ROI_noperm'};
%D.suff = {'t-3000_-2_b-3000_-2500_m_-2500_-1000_Grp_Exp_Subject_orig_cleaned_SPNall_prt_GrpAvCond_gpc_ROI_noperm'};
%D.suff = {'t-5500_1500_b-5500_-5000_m_0_1500_Grp_Exp_Subject_orig_cleaned_trialNmatch_prt_GrpAvCond_gpc_ROI_noperm'};
%D.suff = {'Expectancy\ExpHL\t-3000_-2_b-3000_-2500_m_-2500_-1000_Grp_Exp_Subject_orig_cleaned_SPNall_prt_ExpHL_gpc_ROI'};
%D.suff = {'t-5500_1500_b-5500_-5000_m_0_1500_Grp_Exp_Subject_orig_cleaned_trialNmatch_prt_Exp_gpc_ROI_noperm'};
%D.suff = {'t-3000_-2_b-3000_-2500_m_-2500_-1000_Grp_Exp_Subject_orig_cleaned_SPNall_prt_GrpAvCond_gpc_ROI_noperm',...
    %'t-3000_-2_b-3000_-2500_m_-2500_-1000_Grp_Exp_Subject_orig_cleaned_SPNall_prt_ExpHL_gpc_ROI_noperm'};%,...
    %'Expectancy\ExpHL\t-3000_-2_b-3000_-2500_m_-2500_-1000_Grp_Exp_Subject_orig_cleaned_SPNall_prt_Exp_gpc_ROI'
    %'Group Expectancy\t-3000_-2_b-3000_-2500_m_-2500_-1000_Grp_Exp_Subject_orig_cleaned_SPNall_prt_ExpLowPC_gpc_ROI'
    %'Group Expectancy\t-3000_-2_b-3000_-2500_m_-2500_-1000_Grp_Exp_Subject_orig_cleaned_SPNall_prt_ExpHighPC_gpc_ROI'};
    %'Group Expectancy\t-3000_-2_b-3000_-2500_m_-2500_-1000_Grp_Exp_Subject_orig_cleaned_SPNall_prt_ExpLL_LowPC_gpc_ROI'};
    %'Group Expectancy\t-3000_-2_b-3000_-2500_m_-2500_-1000_Grp_Exp_Subject_orig_cleaned_SPNall_prt_ExpLL_HighPC_gpc_ROI'};
    %''};
%D.suff = {%'t-5500_1500_b-5500_-5000_m_0_1500_Grp_Exp_Subject_orig_cleaned_trialNmatch_prt_GrpAvCond_gpc_ROI_noperm',...
%    't-5500_1500_b-5500_-5000_m_0_1500_Grp_Exp_Subject_orig_cleaned_trialNmatch_prt_Exp_gpc_ROI_noperm'
%    't-5500_1500_b-5500_-5000_m_0_1500_Grp_Exp_Subject_orig_cleaned_trialNmatch_prt_Int_gpc_ROI_noperm'};
    %'Group Expectancy\t-5500_1500_b-5500_-5000_m_0_1500_Grp_Exp_Subject_orig_cleaned_trialNmatch_prt_ExpLowPC_gpc_ROI'};
    %'Group Expectancy\t-5500_1500_b-5500_-5000_m_0_1500_Grp_Exp_Subject_orig_cleaned_trialNmatch_prt_ExpHighPC_gpc_ROI'};
    %''};
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
% input 3: indices of D.suff for dependent variable
% input 4: indices of D.suff for independent variable(s)
% input 5: indices of D.suff for nuisance covariate(s)
%use_bootstrapped=0; 
%regress_images(D,{'wimg','pwimg'},1,2,[],use_bootstrapped)

%% plot weight TOIs
% input 3: indices of D.suff to plot
% input 4: masked time array (for x axis)
% input 5: full time window unmasked
%plot_weight_TOIs(D,{'wimg','pwimg'},1,[-2495:10:-1005],[-3000:2:-2]);
plot_weight_TOIs(D,{'wimg','pwimg'},1,[5:10:1495],[-5500:2:1500]);

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


function create_projections(D,ind) 
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
    try
        fs=PRT.fas.dat(:,:);
    catch % correct path
        [pth nme ext] = fileparts(PRT.fas.dat.fname);
        PRT.fas.dat.fname = fullfile(D.stats_path,D.suff{1},[nme ext]);
        fs=PRT.fas.dat(:,:);
    end
    
    % get indices of feature set
    fsind = PRT.fas.idfeat_img;
    
    % put features into the full array
    fs_full = nan(size(fs,1),length(nonan));
    for ii = 1:size(fs,1)
        fs_full(ii,fsind) = fs(ii,:); 
    end
    
    % reshape fs_full to the dimensions of the weight image
    fsd = reshape(fs_full,size(fs_full,1),size(Wimg,1),size(Wimg,2),size(Wimg,3));
    
    % put feature set in weight image dimensions and convert to 2D
    %fsd = repmat(double(~isnan(Wimg)),1,1,1,size(fs,1));
    %fsd=permute(fsd,[4 1 2 3]);
    %fsd(find(fsd)) = fs(:);
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


function regress_images(D,fields,indDV,indIV,cov,boot)

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
    
    IVind = [indIV,cov];
    
    for i = 1:length(IVind)
        % load IV
        R.IV{i} = D.(fields{f}){IVind(i)};
        V=spm_vol(R.IV{i});
        Viv=spm_read_vols(V(end));
        % vectorise
        iv(:,i) = Viv(:);
    end
    dvi = ~isnan(dv);
    ivi = ~isnan(iv(:,1));
    if sum(dvi)~=sum(ivi)
        %error('dv and iv are different sizes')
        dvi = dvi.*ivi;
        ivi = dvi;
    end
    
    % run regression
    [R.st.stats,rdv] = regress_with_cov(dv(find(dvi)),iv(find(ivi),IVind==indIV),iv(find(ivi),ismember(IVind,cov)));
    
    if boot && strcmp(fields{f},'wimg')
        DVs = dir(fullfile(D.stats_path,[D.pref D.suff{indDV}],'perm_weights*'));
        DV = fullfile(D.stats_path,[D.pref D.suff{indDV}],DVs.name);
        for ii = 1:length(IVind)
            IVs = dir(fullfile(D.stats_path,[D.pref D.suff{IVind(ii)}],'perm_weights*'));
            IV{ii,1} = fullfile(D.stats_path,[D.pref D.suff{IVind(ii)}],IVs.name);
        end
        P = regress_permuted(DV,IV,find(IVind==indIV),find(IVind==cov));
        R.st.permresults = extractfield(P, 'rspacetime'); 
        R.st.pval = sum(abs(R.st.permresults) > abs(R.st.stats(1))) / length(R.st.permresults);
        pID='permuted';
        R.st.CIs = prctile(R.st.permresults,[2.5 97.5]);
    else
        R.st.pval = R.st.stats(3);
        pID='';
        R.st.CIs = [nan,nan];
    end
    figure;scatter(iv(find(ivi),IVind==indIV),rdv)
    title([fields{f} ': space-time r2=' num2str(R.st.stats(1)) ', ' pID ' p=' num2str(R.st.pval) ', CIs = ' num2str(R.st.CIs(1)) ' ' num2str(R.st.CIs(2))])
    
    % regression over space for each time point
    R.s.rtime=nan(size(Vdv,3),1);
    R.s.pval = nan(size(Vdv,3),1);
    R.s.CIs = nan(size(Vdv,3),2);
    ivi=[];
    dvi=[];
    for d = 1:size(Vdv,3)
        ivt=[];
        dvt=[];
        if all(isnan(Vdv(:,:,d)))
            continue;
        end
        dvt = Vdv(:,:,d);
        for i = 1:size(iv,2)
            ivtemp = Viv(:,:,d);
            ivi(:,i) = ~isnan(ivtemp(:));
        end
        dvt=dvt(:);
        dvi = ~isnan(dvt);
        if sum(dvi)~=sum(ivi(1))
            %error('dv and iv are different sizes')
            for i = 1:size(iv,2)
                ivtemp = Viv(:,:,d);
                ivtemp = ivtemp(:);
                ivi(:,i) = ~isnan(ivtemp);
                dvi = dvi.*ivi;
                ivi = dvi;
                ivt(:,i) = ivtemp(find(ivi(:,i)));
            end
        end
        dvt = dvt(find(dvi));
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
            R.s.pval(d) = stats(3);
            pT='';
            R.s.CIs(d,:) = [nan,nan];
        end
    end
    figure;
    rtime = R.s.rtime;
    rtime_corr = nan(size(rtime));
    %if strcmp(pT,'permuted')
        rtime(R.s.pval>0.05)=nan;
        [R.s.pID,R.s.pN] = fdr(R.s.pval,0.05);
        if ~isempty(R.s.pN)
            rtime_corr = rtime;
            rtime_corr(R.s.pval>R.s.pN) = nan;
        end
    %end 
    scatter(1:length(rtime),rtime,'b'); hold on
    scatter(1:length(rtime_corr),rtime_corr,'r'); hold off
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
    for i = 1:length(IVind)
        load(D.prt{IVind(i)});
        ivwt(:,i) = PRT.model.output.weight_ROI{1, 1}(:,end);
        if length(ivwt(:,i))~=length(dvwt)
            error('dv and iv are different sizes')
        end
    end
    
    % run regression
    [R.t.stats,dvwt] = regress_with_cov(dvwt,ivwt(:,IVind==indIV),ivwt(:,ismember(IVind,cov)));
    
    % multiple regression
    %STATS contains:
    %the R-square statistic, 
    %the F statistic 
    %the p value
    %an estimate of the error variance.
    %[b,bint,r,rint,R.t.stats] = regress(dvwt,[ones(length(dvwt),1) ivwt]);
    figure;scatter(ivwt(:,IVind==indIV),dvwt)
    title([fields{f} ': weights over time, r2=' num2str(R.t.stats(1)) ', p=' num2str(R.t.stats(3))])
end
save(fullfile(D.stats_path,[D.pref D.suff{indDV}],'weights_stats.mat'),'R');
end

function P = regress_permuted(DV,IV,indIV,indCOV)
    
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
    % run regression
    stats = regress_with_cov(P(i).dv(P(i).dvi),P(i).iv(P(i).ivi(:,indIV)),P(i).iv(P(i).ivi(:,indCOV)));
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

function [stats,rdv] = regress_with_cov(dv,iv,cov)
% multiple regression, controlling for nuisance covariates
%STATS contains:
%the R-square statistic, 
%the F statistic 
%the p value
%an estimate of the error variance.

% first, get residuals from regressing with covariate
r=[];
if ~isempty(cov)
    for c = 1:size(cov,2)
        if c==1;
            rdv = dv;
        else
            rdv = r;
        end
        [~,~,r,~,~] = regress(rdv,[ones(size(cov,1),1) cov(:,c)]);
    end
end

% then, regress residuals with IV of interest
if ~isempty(r);
    rdv = r;
else
    rdv = dv;
end
[b,bint,r,rint,stats] = regress(rdv,[ones(size(iv,1),1) iv]);
end

function plot_weight_TOIs(D,field,ind,xval,tp)

DV = D.(field{1}){ind};
pth=fileparts(DV);
load(D.prt{ind});
dvwt = PRT.model.output.weight_ROI{1, 1}(:,end);

% for plotting colours in variance range
%figure;
%set(gcf, 'Position', [100, 100, 1000, 300])
%scatter(xval,dvwt,[],[0.5 0.5 0.5],'filled'); hold on
%sds = [0 1 2];
%cols = {[0 0 1],[0 0.8 0],[1 0 0]};
%cols = linspace(1,10,length(sds));


%plot([min(xval),max(xval)],[thresh,thresh],'r')
%title([field ': weights over time'])
%xlabel('Time (ms)')
%ylabel('Weight')

si_peak = find(dvwt==max(dvwt));
xval_peak = xval(si_peak);
disp(['peak value: ' num2str(xval_peak)]);

C=zeros(1,length(xval));
sds = [-3:1:3];
for s = 1:length(sds)
    sd=sds(s);
    thresh = mean(dvwt)+sd*std(dvwt);
    si = find(dvwt>thresh);
    C(si) = sd;
end
figure
set(gcf, 'Position', [100, 100, 1000, 300])
x = xval;
y = [0.5];
%C = 1:1500;
imagesc(x,y,C)
colormap(jet)
colorbar

% get weight matrix imge 
for f = 1:length(field)
    W = D.(field{f}){ind};
    if ~exist(W)
        disp('no projections');
        continue
    end
    V=spm_vol(W);
    R(f).Vdv=spm_read_vols(V(end));
    
    % plot temporal summaries
    
    Vdv = reshape(R(f).Vdv,size(R(f).Vdv,1)*size(R(f).Vdv,2),size(R(f).Vdv,3));
    incr = (xval(2)-xval(1))/2; % increment in TOI image
    for x = 1:length(xval)
        xv = xval(x)-incr+1:xval(x)+incr;
        tw=find(tp>=xv(1) & tp<=xv(end));
        temp=Vdv(:,tw);
        R(f).TOImean(x,1) = nanmean(abs(temp(:)));
        R(f).TOIstd(x,1) = nanstd(temp(:));
    end
end

f1=figure
set(gcf, 'Position', [100, 100, 1000, 300])
s=scatter(xval,dvwt*100/sum(dvwt),100,C,'filled'); 
s.MarkerEdgeColor = [0 0 0];
s.LineWidth=0.5;
colormap(jet)
caxis([min(sds),max(sds)])
%title([field ': weights over time'])
xlabel('Time (ms)')
ylabel('Weight (% contribution)')
%axis off
tightfig(f1)
colorbar

% get weight matrix img
T=table;
for f = 1:length(field)
    
    W = D.(field{f}){ind};
    if ~exist(W)
        disp('no projections');
        continue
    end
    
    dat = R(f).TOImean;
    si_peak = find(dat==max(dat));
    xval_peak = xval(si_peak);
    for s = 1:length(sds)
        sd=sds(s);
        thresh = mean(dat)+sd*std(dat);
        si = find(dat>thresh);
        C(si) = sd;
    end
    fig(f)=figure
    set(gcf, 'Position', [100, 100, 1000, 300])
    s=scatter(xval,dat*100/sum(dat),100,C,'filled'); 
    s.MarkerEdgeColor = [0 0 0];
    s.LineWidth=0.5;
    colormap(jet)
    caxis([min(sds),max(sds)])
    title([field{f} ': weights over time - means'])
    xlabel('Time (ms)')
    ylabel('Weight (% contribution)')
    %axis off
    %tightfig(fig(f))
    colorbar
    
    % Tabulate
    T.time = xval';
    T.([field{f} '_val']) = dat;
    T.([field{f} '_stdev']) = C';
    
    % plot max TOI topo/projection
    incr = (xval(2)-xval(1))/2; % increment in TOI image
    tw=find(tp>xval_peak-incr & tp<xval_peak+incr);
    for f2 = 1:length(field)
        if length(R)>=f2
            topo = nanmean(R(f2).Vdv(:,:,tw),3);
            figure
            set(gcf, 'Position', [400, 400, 350, 300])
            %imagesc(rot90(topo))
            pcolor(rot90(topo,3)), shading interp
            axis off
            colormap(jet)
            title([field{f} ': time of max ' field{f2} ' for image ' num2str(xval_peak) 'ms'])
            colorbar
        end
    end
end

writetable(T,fullfile(pth,'weight_timeseries.xlsx'));
cd(pth)

[sdv,si] = sort(dvwt,'descend');
for s = 1:length(sds)
    sd1=sds(s);
    sd2=sds(min(s+1,length(sds)));
    thresh1 = mean(dvwt)+sd1*std(dvwt);
    if sd2>sd1
        thresh2 = mean(dvwt)+sd2*std(dvwt);
    else
        thresh2 = Inf;
    end
    si = find(dvwt>thresh1 & dvwt<thresh2);
    xval_max{s} = nan(size(xval));
    dvwt_max{s} = nan(size(dvwt));
    xval_max{s}(si) = xval(si);
    dvwt_max{s}(si) = dvwt(si);
%    scatter(xval_max{s},dvwt_max{s},[],cols{s},'filled')
    disp(['max values (SD: ' num2str(sd1) '): ' num2str(xval_max{s}(~isnan(xval_max{s})))]);
end

end