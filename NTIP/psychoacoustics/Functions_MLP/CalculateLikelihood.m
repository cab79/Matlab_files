function ll=CalculateLikelihood(x, responses, alpha, beta, gamma, lambda)
%
% ll = CalculateLikelihood(x, responses, alpha, beta, gamma, lambda)
%
% This function returns the likelihood of a given psychometric function.
%
% - X: an array containing the stimuli levels presented thus far;
% - RESPONSES: an array containing the subject's responses collected thus
% far. Responses must be coded as "0" (wrong or "no"), % "1" (correct, or
% "yes");
% - ALPHA: the midpoint of the psychometric function;
% - BETA: the slope of the psychometric function;
% - GAMMA: the lower limit of the psychometric function expressed as a
% proportion;
% - LAMBDA: the subject's lapse rate expressed as a proportion.

warning off

ll = 0;
for i=1:length(x)
    p=Logistic(x(i), alpha, beta, gamma, lambda);
    if responses(i)==1
        ll=ll+log(p);
    else
        ll=ll+log(1-p);
    end;
end;
