function out = bayesreg_crossval(X,y,S,groupvec)

if S.zscore
    [X, X_mean, X_stand_de] = zscore(X, [], 1);
    [y, y_mean, Y_stand_de] = zscore(y);
end

N=size(X,1);
K=S.brr.folds;
if K
    x_val_in  = crossvalind('Kfold', N, K);
else
    % if no folds, run all trials for both training and testing
    x_val_in = ones(N,1);
    K=1;
end
for i = K : -1 : 1
    
    val_in    = x_val_in == i;
    es_in     = ~val_in;
    if ~any(es_in); es_in=val_in; end
    
    if S.brr.usegroups
        grps = unique(groupvec);
        if length(grps)>1
            for g = grps
                groups{g} = find(groupvec==g);
            end
            [beta, beta0, stt(i)] = bayesreg(X(es_in,:),y(es_in),S.brr.model,S.brr.prior,'nsamples',S.brr.nsamples,'burnin',S.brr.burnin,'thin',S.brr.thin,'display',false,'waic', S.brr.waic, 'groups',groups);
        else
            [beta, beta0, stt(i)] = bayesreg(X(es_in,:),y(es_in),S.brr.model,S.brr.prior,'nsamples',S.brr.nsamples,'burnin',S.brr.burnin,'thin',S.brr.thin,'display',false,'waic', S.brr.waic);
        end
    else
        [beta, beta0, stt(i)] = bayesreg(X(es_in,:),y(es_in),S.brr.model,S.brr.prior,'nsamples',S.brr.nsamples,'burnin',S.brr.burnin,'thin',S.brr.thin,'display',false,'waic', S.brr.waic);
    end
    [pred, predstt(i)] = br_predict(X(val_in,:), beta, beta0, stt(i), 'ytest', y(val_in), 'CI', [2.5, 97.5], 'display', false);
    
    logl(i) = stt(i).modelstats.logl; 
    waic(i) = stt(i).modelstats.waic; 
    r2(i) = stt(i).modelstats.r2; 
end

out.muB = mean([stt(:).muB],2);
out.muSigma2 = mean([stt(:).muSigma2]);
out.waic = mean(waic);
out.logl = mean(logl);
out.r2 = mean(r2);
out.neglike = mean([predstt(:).neglike]);
out.r2test = mean([predstt(:).r2]);
