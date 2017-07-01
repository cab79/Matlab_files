function numcomp = numcompeig(EEG);

% determine numcomponent by doing an eig on the covariance matrix
covar = zeros(EEG.nbchan);
for itrial = 1:EEG.trials
    currtrial = EEG.data(:,:,itrial);
    covar = covar + currtrial*currtrial.';
end
[V, D] = eig(covar);
D = sort(diag(D),'descend');
D = D ./ sum(D);
Dcum = cumsum(D);
numcomp = find(Dcum>.99,1,'first');