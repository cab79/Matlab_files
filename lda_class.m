function out = lda_class(dat,groups,nfold,ndec,perm)
% multivariate classification between two groups
% using linear discriminant analysis

dat=round(dat,ndec);
rm=all(isnan(dat),1) | all(diff(dat)==0); % remove nans and constants
dat(:,rm)=[];

if isempty(dat)
    out.empty=1;
    return
end

try
    lda = fitcdiscr(dat,groups,'DiscrimType','linear');
catch
    disp('trying pseudolinear LDA')
    lda = fitcdiscr(dat,groups,'DiscrimType','pseudolinear');
end
ldaClass = resubPredict(lda);
% misclassification error (the proportion of misclassified observations) on the training set.
ldaResubErr = resubLoss(lda);
% confusion matrix
[ldaResubCM,grpOrder] = confusionmat(groups,ldaClass);
%Estimate the true test error for LDA using 10-fold stratified cross-validation.
cp = cvpartition(groups,'KFold',nfold);
cvlda = crossval(lda,'CVPartition',cp);
out.ldaCVErr = kfoldLoss(cvlda);

if perm
    for i = 1:perm
        disp(['LDA permutation ' num2str(i) '/' num2str(perm)])
        groups = groups(randperm(length(groups)));
        try
            lda = fitcdiscr(dat,groups,'DiscrimType','linear');
        catch
            disp('trying pseudolinear LDA')
            lda = fitcdiscr(dat,groups,'DiscrimType','pseudolinear');
        end
        ldaClass = resubPredict(lda);
        % misclassification error (the proportion of misclassified observations) on the training set.
        ldaResubErr = resubLoss(lda);
        % confusion matrix
        [ldaResubCM,grpOrder] = confusionmat(groups,ldaClass);
        %Estimate the true test error for LDA using 10-fold stratified cross-validation.
        cp = cvpartition(groups,'KFold',nfold);
        cvlda = crossval(lda,'CVPartition',cp);
        out.perm.ldaCVErr(i) = kfoldLoss(cvlda);
    end
    out.pval = sum(out.perm.ldaCVErr<out.ldaCVErr)/perm;
    
end

% clear 'ClassificationDiscriminant' ...
%     'classreg.learning.modelparams.DiscriminantParams' ...
%     'classreg.learning.modelparams.EnsembleParams' ...
%     'classreg.learning.generator.Partitioner' ...
%     'classreg.learning.classif.CompactClassificationDiscriminant' ...
%     'classreg.learning.FitTemplate' ...
%     'update' ...
%     'cvpartition' ...
%     'classreg.learning.internal.ClassLabel' ...
%     'classreg.learning.modifier.BlankModifier' ...
%     'classreg.learning.combiner.WeightedAverage' ...
%     'classreg.learning.impl.CompactEnsembleImpl'
%     
% [M,X,C] = inmem;
% C