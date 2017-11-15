function DCM_PEB(GCMfile,DCMdir,run_num,fullmodel,loadorsave,excl_full)
dbstop if error
load(fullfile(DCMdir,GCMfile));
cd(DCMdir)
DCM = GCM{1,fullmodel};
% The following section contains the key analyses
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

% Bayesian model reduction (avoiding local minima over models)
%==========================================================================
% spm_dcm_bmr operates on different DCMs of the same data (rows) to find
% the best model. It assumes the full model - whose free-parameters are
% the union (superset) of all free parameters in each model - has been
% inverted, and from that model creates reduced versions for the other
% models. If the other models have already been inverted, these are ignored
% and only the full model inversion is used to derive the others.
% GCM is a {Nsub x Nmodel} cell array of DCM filenames or model structures  
%         of Nsub subjects, where each model is reduced independently 
rname=['RCM_fit' num2str(run_num) '.mat'];
if loadorsave && exist(fullfile(DCMdir,rname),'file')
    load(fullfile(DCMdir,rname));
else
    RCM   = spm_dcm_bmr(GCM);
    if iscell(RCM{1,1}); RCM = vertcat(RCM{:});end;
    save(rname,'RCM','-v7.3');
end

% remove the full model, assuming it is biologically implausible
if excl_full
    GCM(:,fullmodel) =[];
    RCM(:,fullmodel) =[];
end

% hierarchical (empirical Bayes) model reduction:
% optimises the empirical priors over the parameters of a set of first level DCMs, using second level or
% between subject constraints specified in the design matrix X.
% See: https://en.wikibooks.org/wiki/SPM/Parametric_Empirical_Bayes_(PEB)
%==========================================================================
%[PEB,P]   = spm_dcm_peb(P,M,field)
% 'Field' refers to the parameter fields in DCM{i}.Ep to optimise [default: {'A','B'}]
% 'All' will invoke all fields. This argument effectively allows 
%  one to specify the parameters that constitute random effects. 
% If P is an an {N [x M]} array, will return one PEB for all
% models/subjects
% If P is an an (N x M} array, will return a number of PEBs as a 
% {1 x M} cell array.
pname=['PCM_peb' num2str(run_num) '.mat'];
if loadorsave && exist(fullfile(DCMdir,pname),'file')
    load(fullfile(DCMdir,pname));
else
    [peb,PCM] = spm_dcm_peb(RCM,[],{'A','B'});
    save(pname,'PCM','peb','-v7.3');
end

% BMA - first level
%--------------------------------------------------------------------------
% Bayesian model averages, weighting each model by its marginal likelihood pooled over subjects
try bma  = spm_dcm_bma(GCM); end;
rma  = spm_dcm_bma(RCM);
pma  = spm_dcm_bma(PCM);

% BMC/BMA - second level
%==========================================================================

% second level model
%--------------------------------------------------------------------------
M    = struct('X',Xb);

% BMC - search over first and second level effects
%--------------------------------------------------------------------------
[BMC,PEB] = spm_dcm_bmc_peb(RCM,M,{'A','B'});

% BMA - exhaustive search over second level parameters
%--------------------------------------------------------------------------
PEB.gamma = 1/128;
BMA       = spm_dcm_peb_bmc(PEB);


% posterior predictive density and LOO cross validation
%==========================================================================
if length(unique(Xb))>1
    spm_dcm_loo(RCM(:,1),Xb,{'B'});
end

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



% extract and plot results
%==========================================================================
[Ns Nm] = size(GCM);
clear P Q R
for i = 1:Ns
    
    % data - over subjects
    %----------------------------------------------------------------------
    for nt = 1:length(GCM{i,1}.xY.y)
        Y(:,i,nt) = GCM{i,1}.xY.y{nt}*DCM.M.U(:,1);
    end
    pst = GCM{i,1}.xY.pst;
    
    % Parameter averages (over models)
    %----------------------------------------------------------------------
    try Q(:,i,1) = full(spm_vec(bma.SUB(i).Ep));end;
    Q(:,i,2) = full(spm_vec(rma.SUB(i).Ep));
    Q(:,i,3) = full(spm_vec(pma.SUB(i).Ep));
    
    % Free energies
    %----------------------------------------------------------------------
    for j = 1:Nm
        try F(i,j,1) = GCM{i,j}.F - GCM{i,1}.F;end;
        F(i,j,2) = RCM{i,j}.F - RCM{i,1}.F;
        F(i,j,3) = PCM{i,j}.F - PCM{i,1}.F;
    end
    
end


% classical inference
%==========================================================================

% indices to plot parameters
%--------------------------------------------------------------------------
pC    = DCM.M.pC;
iA    = spm_find_pC(pC,pC,'A');
iB    = spm_find_pC(pC,pC,'B');


% classical inference of second level
%--------------------------------------------------------------------------
i   = [iA;iB];
%CVA = spm_cva(Q(i,:,1)',Xb,[],[0 1]'); 
%CVA = spm_cva(Q(i,:,2)',Xb,[],[0 1]'); 
CVA = spm_cva(Q(i,:,3)',Xb,[],[0 1]'); 

spm_figure('GetWin','CVA');clf
bar(spm_en([Xb(:,2) CVA.v]))
title('Canonical variate','FontSize',16)
xlabel('parameter'), ylabel('weight'), axis square
legend({'RFX','True'})

% plot data
%==========================================================================
spm_figure('GetWin','Figure 1');clf

p =  (Xb(:,2) > 0);
q = ~(Xb(:,2) > 0);

subplot(2,2,1)
if any(q)
    col = {'r','r:','r--'};
    for np = 1:size(Y,3)
        plot(pst,Y(:,q,np),col{np}),  hold on
    end
end
if any(p)
    col = {'b','b:','b--'};
    for np = 1:size(Y,3)
        plot(pst,Y(:,p,np),col{np}),  hold on
    end
end
hold off
xlabel('pst'), ylabel('response'), title('Group data','FontSize',16)
axis square, spm_axis tight

try
    subplot(2,2,2)
    if any(q); plot(pst,Y(:,q,1) - Y(:,q,2),'r'), hold on; end
    if any(p); plot(pst,Y(:,p,1) - Y(:,p,2),'b'), hold off; end
    xlabel('pst'), ylabel('differential response'), title('Difference waveforms','FontSize',16)
    axis square, spm_axis tight
end

i = spm_fieldindices(DCM.Ep,'B{1}(1,1)');
j = spm_fieldindices(DCM.Ep,'B{1}(2,2)');

subplot(2,2,3)
plot(Q(i,q,3),Q(j,q,3),'.r','MarkerSize',24), hold on
plot(Q(i,p,3),Q(j,p,3),'.b','MarkerSize',24), hold off
xlabel('B{1}(1,1)'), ylabel('B{1}(2,2)'), title('Group effects PMA','FontSize',16)
axis square

i = spm_fieldindices(DCM.Ep,'B{1}(3,3)');
j = spm_fieldindices(DCM.Ep,'B{1}(4,4)');

subplot(2,2,4)
plot(Q(i,q,3),Q(j,q,3),'.r','MarkerSize',24), hold on
plot(Q(i,p,3),Q(j,p,3),'.b','MarkerSize',24), hold off
xlabel('B{1}(3,3)'), ylabel('B{1}(4,4)'), title('Group effects PMA','FontSize',16)
axis square
print('fig1','-dpng') 

% plot results: Bayesian model reduction vs. reduced models
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2: Pooled free energy on Right'); clf

occ = 512;
try
    f   = F(:,:,1); f = f - max(f(:)) + occ; f(f < 0) = 0;
    subplot(3,2,1), imagesc(f)
    xlabel('model'), ylabel('subject'), title('Free energy (FFX)','FontSize',16)
    axis square

    f   = sum(f,1); f  = f - max(f) + occ; f(f < 0) = 0;
    subplot(3,2,3), bar(f), xlabel('model'), ylabel('Free energy'), title('Free energy (FFX)','FontSize',16)
    spm_axis tight, axis square

    p   = softmax(f'); [m,i] = max(p);
    subplot(3,2,5), bar(p)
    text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','k','FontSize',8)
    xlabel('model'), ylabel('probability'), title('Posterior (FFX)','FontSize',16)
    axis([0 (length(p) + 1) 0 1]), axis square
end

occ = 128;
f   = F(:,:,2); f = f - max(f(:)) + occ; f(f < 0) = 0;
subplot(3,2,2), imagesc(f)
xlabel('model'), ylabel('subject'), title('Free energy (BMR)','FontSize',16)
axis square

f   = sum(f,1); f  = f - max(f) + occ; f(f < 0) = 0;
subplot(3,2,4), bar(f), xlabel('model'), ylabel('Free energy'), title('Free energy (BMR)','FontSize',16)
spm_axis tight, axis square

p   = softmax(f'); [m,i] = max(p);
subplot(3,2,6), bar(p)
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','k','FontSize',8)
xlabel('model'), ylabel('probability'), title('Posterior (BMR)','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square

print('fig2','-dpng') 

% a more detailed analysis of Bayesian model comparison for BMR
% compares subjects with best and worst evidence for the full model
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3'); clf

f   = F(:,:,2); f = f - max(f(:)) + occ; f(f < 0) = 0;
subplot(3,2,1), imagesc(f)
xlabel('model'), ylabel('subject'), title('Free energy (BMR)','FontSize',16)
axis square

pp  = softmax(f')';
subplot(3,2,2), imagesc(pp)
xlabel('model'), ylabel('subject'), title('Model posterior (BMR)','FontSize',16)
axis square

[p,i] = max(pp(:,fullmodel));
[p,j] = min(pp(:,fullmodel));
stri  = sprintf('Subject %i',i);
strj  = sprintf('Subject %i',j);

subplot(3,2,3), bar(pp(i,:))
xlabel('model'), ylabel('probability'), title(stri,'FontSize',16)
axis square, spm_axis tight

subplot(3,2,4), bar(pp(j,:))
xlabel('model'), ylabel('probability'), title(strj,'FontSize',16)
axis square, spm_axis tight

k   = spm_fieldindices(DCM.Ep,'B');
qE  = RCM{i,1}.Ep.B{1}; qE = spm_vec(qE);
qC  = RCM{i,1}.Cp(k,k); qC = diag(qC);
qE  = qE(find(qC));
qC  = qC(find(qC));

subplot(3,2,5), spm_plot_ci(qE,qC),
xlabel('parameter (B) for each connection'), ylabel('expectation'), title('Parameters','FontSize',16)
axis square, a = axis;

qE  = RCM{j,1}.Ep.B{1}; qE = spm_vec(qE);
qC  = RCM{j,1}.Cp(k,k); qC = diag(qC);
qE  = qE(find(qC));
qC  = qC(find(qC));

subplot(3,2,6), spm_plot_ci(qE,qC), 
xlabel('parameter (B)'), ylabel('expectation'), title('Parameters','FontSize',16)
axis square, axis(a);

print('fig3','-dpng') 

% a more detailed analysis of Bayesian model comparison for PEB
% compares subjects with best and worst evidence for the full model
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 4'); clf

f   = F(:,:,3); f = f - max(f(:)) + occ; f(f < 0) = 0;
subplot(3,2,1), imagesc(f)
xlabel('model'), ylabel('subject'), title('Free energy (PCM)','FontSize',16)
axis square

pp  = softmax(f')';
subplot(3,2,2), imagesc(pp)
xlabel('model'), ylabel('subject'), title('Model posterior (PCM)','FontSize',16)
axis square

[p,i] = max(pp(:,fullmodel));
[p,j] = min(pp(:,fullmodel));
stri  = sprintf('Subject %i',i);
strj  = sprintf('Subject %i',j);

subplot(3,2,3), bar(pp(i,:))
xlabel('model'), ylabel('probability'), title(stri,'FontSize',16)
axis square, spm_axis tight

subplot(3,2,4), bar(pp(j,:))
xlabel('model'), ylabel('probability'), title(strj,'FontSize',16)
axis square, spm_axis tight

k   = spm_fieldindices(DCM.Ep,'B');
qE  = PCM{i,1}.Ep.B{1}; qE = spm_vec(qE);
qC  = PCM{i,1}.Cp(k,k); qC = diag(qC);
qE  = qE(find(qC));
qC  = qC(find(qC));

subplot(3,2,5), spm_plot_ci(qE,qC),
xlabel('parameter (B) for each connection'), ylabel('expectation'), title('Parameters','FontSize',16)
axis square, a = axis;

qE  = PCM{j,1}.Ep.B{1}; qE = spm_vec(qE);
qC  = PCM{j,1}.Cp(k,k); qC = diag(qC);
qE  = qE(find(qC));
qC  = qC(find(qC));

subplot(3,2,6), spm_plot_ci(qE,qC), 
xlabel('parameter (B)'), ylabel('expectation'), title('Parameters','FontSize',16)
axis square, axis(a);

print('fig4','-dpng') 

PCM{i,1}.K = DCM.K;
PCM{i,1}.H = DCM.H;
PCM{i,1}.R = DCM.R;
PCM{j,1}.K = DCM.K;
PCM{j,1}.H = DCM.H;
PCM{j,1}.R = DCM.R;
bs=spm_figure('GetWin','Modes - best subject'); %clf
spm_dcm_erp_results(PCM{i,1},'ERPs (mode)',bs);
ws=spm_figure('GetWin','Modes - worst subject'); %clf
spm_dcm_erp_results(PCM{j,1},'ERPs (mode)',ws);


% first level parameter estimates and Bayesian model averages
% correlates parameters of the BMA across the 3 schemes (cf "Correlations"
% figure which is on the true model, not the average across models)
%--------------------------------------------------------------------------
spm_figure('GetWin','Summed free energy');clf, ALim = 1/2;
try 
    p   = spm_softmax(sum(F(:,:,1))');
    subplot(3,1,1), bar(p),[m,i] = max(p);
    text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','k','FontSize',8)
    xlabel('model'), ylabel('probability'), title('Posterior (FFX)','FontSize',16)
    axis([0 (length(p) + 1) 0 1]), axis square
end

p   = spm_softmax(sum(F(:,:,2))');
subplot(3,1,2), bar(p),[m,i] = max(p);
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','k','FontSize',8)
xlabel('model'), ylabel('probability'), title('Posterior (BMR)','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square

p   = spm_softmax(sum(F(:,:,3))');
subplot(3,1,3), bar(p),[m,i] = max(p);
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','k','FontSize',8)
xlabel('model'), ylabel('probability'), title('Posterior (PEB)','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square

print('fig5','-dpng') 

% random effects Bayesian model comparison
%==========================================================================
spm_figure('GetWin','Figure 6: random effects Bayesian model comparison');clf
[~,~,xp] = spm_dcm_bmc(RCM);

p   = BMC.Pw;
subplot(2,2,1), bar(p),[m,i] = max(p);
text(i - 1/4,m/2,sprintf('model %d, %-2.0f%%',i,m*100),'Color','k','FontSize',8)
xlabel('model'), ylabel('RCM posterior probability'), title('Random parameter effects','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square
save BMC_Pw p

p   = xp;
subplot(2,2,2), bar(p),[m,i] = max(p);
text(i - 1/4,m/2,sprintf('model %d, %-2.0f%%',i,m*100),'Color','k','FontSize',8)
xlabel('model'), ylabel('RCM exceedance probability'), title('Random model effects','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square
save RCM_xp p

print('fig6','-dpng') 

return



% Notes
%==========================================================================
hE    = linspace(-4,4,16);
hC    = 1;
clear Eh HF
for i = 1:length(hE)
    M.X     = X(:,1:2);
    M.hE    = hE(i);
    M.hC    = 1;
    PEB     = spm_dcm_peb(GCM(:,1),M,{'A','B'});
    HF(i)   = PEB.F;
    Eh(:,i) = PEB.Eh;
    
end

subplot(2,2,1)
plot(hE,HF - max(HF))
title('Free-energy','FontSize',16)
xlabel('Prior')
subplot(2,2,2)
plot(hE,Eh)
xlabel('Prior')
title('log-preciion','FontSize',16)

