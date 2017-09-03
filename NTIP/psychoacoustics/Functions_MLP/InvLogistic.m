function x = InvLogistic(p_target, alpha, beta, gamma, lambda)
%
% x = InvLogistic(p_target, alpha, beta, gamma, lambda)
%
% The inverse of the logistic psychometric function
%
% - P_TARGET: the performance we are tracking expressed as a proportion;
% - ALPHA: the midpoint of the psychometric function. It is the value
% of the stimulus where the function crosses the mean probability between
% gamma and lambda (i.e., .75 for a 2AFC, or .66 for a 3AFC, etc.);
% - BETA: the slope of the function. It is the rate of change in the
% subject's performance with stimulus level;
% - GAMMA: the lower limit of the psychometric function expressed as a
% proportion. It ranges from 0 to 1. The meaning of gamma is different for
% yes/no and nAFC tasks. In nAFC tasks gamma corresponds to chance level
% (i.e., 0.5 in a 2AFC task). In yes/no tasks gamma corresponds to the
% subject's false alarm rate;
% - LAMBDA: the subject's lapse rate expressed as a proportion.

if nargin<5
    lambda=0;
end;

x=alpha-(1/beta)*log(((1-lambda-gamma)/(p_target-gamma))-1);