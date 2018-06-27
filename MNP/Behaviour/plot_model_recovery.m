function MNP_plot_model_recovery(D,S)

close all
clear all
cd('C:\Data\MNP\Pilots\VLTv3\processed')

basename = 'evidence_model';

% for each simulated model
sim_mod = 3:7;
for m = 1:length(sim_mod)
    load([basename num2str(sim_mod(m))]);
    
    % mean evidence for each fitted model over repetitions
    LME(:,m) = mean(evi.LME,2); % rows: fitted models; columns: simulated models
    AIC(:,m) = mean(evi.AIC,2);
    BIC(:,m) = mean(evi.BIC,2);
end


figure;imagesc(LME), title('LME'), xlabel('simulated'), ylabel('fitted')
set(gca, 'XTick', 1:length(sim_mod), 'XTickLabel', sim_mod) 
set(gca, 'YTick', 1:length(sim_mod), 'YTickLabel', sim_mod) 
figure;imagesc(AIC), title('AIC'), xlabel('simulated'), ylabel('fitted')
set(gca, 'XTick', 1:length(sim_mod), 'XTickLabel', sim_mod) 
set(gca, 'YTick', 1:length(sim_mod), 'YTickLabel', sim_mod) 
figure;imagesc(BIC), title('BIC'), xlabel('simulated'), ylabel('fitted')
set(gca, 'XTick', 1:length(sim_mod), 'XTickLabel', sim_mod) 
set(gca, 'YTick', 1:length(sim_mod), 'YTickLabel', sim_mod) 