function p = Logistic(x, alpha, beta, gamma, lambda)
%
% p = Logistic(x, alpha, beta, gamma, lambda)
%
% The logistic psychometric function
%
% - X:  level of the stimulus;
% - ALPHA: midpoint of the psychometric function. It is the value
% of the stimulus where the function crosses the mean probability between
% gamma and lambda (i.e., .75 for a 2AFC, or .66 for a 3AFC, etc.);
% - BETA: slope of the function. It is the rate of change in the
% subject's performance with stimulus level;
% - GAMMA: lower limit of the psychometric function expressed as a
% proportion. It ranges from 0 to 1. The meaning of gamma is different for
% yes/no and nAFC tasks. In nAFC tasks gamma corresponds to chance level
% (i.e., 0.5 in a 2AFC task). In yes/no tasks gamma corresponds to the
% false alarm rate of the subject;
% - LAMBDA: the subject's lapse rate expressed as a proportion.

if nargin<5
    lambda=0;
end;

p=gamma+((1-lambda-gamma)*(1/(1+exp(beta*(alpha-x)))));