function [out,S,X] = bayesreg_crossval(X,Y,S,groupvec)
% Y can be a matrix with multiple columns, each is analysed separately.

for s = 1:size(Y,2)
    y = Y(:,s);

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

        if (any(isnan(y(es_in))))
            stt(i).muB = NaN*ones(size(X,2),1); 
            stt(i).muSigma2 = NaN; 
            logl(i) = NaN; 
            waic(i) = NaN; 
            r2(i) = NaN; 
            predstt(i).neglike = NaN; 
            predstt(i).r2 = NaN; 
        else
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
    end

    out(s).muB = nanmean([stt(:).muB],2);
    out(s).muSigma2 = nanmean([stt(:).muSigma2]);
    out(s).waic = nanmean(waic);
    out(s).logl = nanmean(logl);
    out(s).r2 = nanmean(r2);
    out(s).neglike = nanmean([predstt(:).neglike]);
    out(s).r2test = nanmean([predstt(:).r2]);
    out(s).xmean = X_mean;
    out(s).xstd = X_stand_de;
    out(s).ymean = y_mean;
    out(s).ystd = Y_stand_de;
end