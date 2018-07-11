function HGF_plot_model_recovery(S)

cd(S.path.prep)

basename = 'evidence_model';

% for each simulated model
sim_mod = S.plot_sim_models;
for m = 1:length(sim_mod)
    files = dir([basename num2str(sim_mod(m)) '*']);
    if length(files)~=1
        error('Wrong number of evidence output files in folder')
    end
    load(files.name);
    
    % mean evidence for each fitted model over repetitions
    LME(:,m) = mean(evi.LME,2); % rows: fitted models; columns: simulated models
    AIC(:,m) = mean(evi.AIC,2);
    BIC(:,m) = mean(evi.BIC,2);
end

figure;imagesc(LME), title('LME'), xlabel('simulated'), ylabel('fitted')
set(gca, 'XTick', 1:length(sim_mod), 'XTickLabel', sim_mod) 
set(gca, 'YTick', 1:length(S.plot_fit_models), 'YTickLabel', S.plot_fit_models) 
figure;imagesc(AIC), title('AIC'), xlabel('simulated'), ylabel('fitted')
set(gca, 'XTick', 1:length(sim_mod), 'XTickLabel', sim_mod) 
set(gca, 'YTick', 1:length(S.plot_fit_models), 'YTickLabel', S.plot_fit_models) 
figure;imagesc(BIC), title('BIC'), xlabel('simulated'), ylabel('fitted')
set(gca, 'XTick', 1:length(sim_mod), 'XTickLabel', sim_mod) 
set(gca, 'YTick', 1:length(S.plot_fit_models), 'YTickLabel', S.plot_fit_models) 