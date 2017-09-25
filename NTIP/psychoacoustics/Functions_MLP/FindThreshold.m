function [level, FA] = FindThreshold(p_target, x, responses, alpha, beta, gamma, lambda)
%
% [level, FA] = FindThreshold (p_target, x, responses, alpha, beta, gamma, lambda)
%
% This function first looks for the most likely psychometric function
% within the range "firstpsyfun", "lastpsyfun" of the vector alpha.
% Successively, it calculates the stimulus level for the next trial at the
% desired level of performance p_target, and returns this value to the user
% together with an estimate of the subject's false alarm rate FA. If the
% task is a nAFC, FA is equal to chance level. 
%
% - P_TARGET: the subject's performance we are targeting; it is the point
% of the psychometric function we are tracking; 
% - X: an array containing the stimuli levels presented thus far;
% - RESPONSES: an array containing the subject's responses collected thus
% far. Responses must be coded as "0" (wrong or "no"), "1" (correct, or
% "yes");
% - ALPHA: an array containing the midpoins of the psychometric functions
% that will drive the threshold search of the maximum likelihood procedure.
% - BETA: the slope of the psychometric function.
% - GAMMA: the lower limit of the psychometric function expressed as a
% proportion. In nAFC tasks GAMMA is a single value (e.g., 0.5 for 2AFC).
% In yes/no tasks it can be an array that contains values of false alarm
% rates.
% - LAMBDA: the subject's lapse rate expressed as a proportion.

ll=zeros(length(alpha), length(gamma));

% calculate the likelihood of each psychometric function
for i=1:length(alpha)
    for j=1:length(gamma)
        ll(i, j)=CalculateLikelihood(x, responses, alpha(i), beta, gamma(j), lambda);
    end;
end;

% find the most likely psychometric function
[i, j]=find(ll==max(max(ll)));
if length(i)+length(j) > 2
    i = i(1);
    j = j(1);
end;

% calculate the level of the stimulus at p_target performance
level=InvLogistic(p_target, alpha(i), beta, gamma(j), lambda);
FA=gamma(j);