function [h,p,varargout] = VBA_groupBMC_btwGroups_CAB(Ls,options)
% test for between-groups difference in model frequencies
% function [h,p] = VBA_groupBMC_btwGroups(Ls,options)
% IN:
%   - Ls: {nmXns_1, nmXns_2} array of log-model evidences matrices of each group (nm models; ns_g subjects in the group).
%   - options: a structure containing the following fields:
%       .DisplayWin: flag for display window
%       .verbose: flag for summary statistics display
%       .families: a cell array of size nf, which contains the indices of
%       the models that belong to each of the nf families.
% OUT:
%   - h: test decision about a difference between the group (rejection of
%        the null hypothesis of equality)
%   - p: the posterior probability that the two groups have the same model
%        frequencies

if nargin < 2
    options = {};
end
% one frequency for all
if length(Ls)==2
    L = [Ls{1} Ls{2}];
elseif length(Ls)==3
    L = [Ls{1} Ls{2} Ls{3}];
end
[posterior,out] = VBA_groupBMC(L,options);
Fe = out.F(end);
varargout = {out};

% separate frequencies
if length(Ls)==2
    [posterior1,out1] = VBA_groupBMC(Ls{1},options);
    [posterior2,out2] = VBA_groupBMC(Ls{2},options);
    Fd = out1.F(end) + out2.F(end);
    varargout = {[varargout,{out1},{out2}]};
elseif length(Ls)==3
    [posterior1,out1] = VBA_groupBMC(Ls{1},options);
    [posterior2,out2] = VBA_groupBMC(Ls{2},options);
    [posterior3,out3] = VBA_groupBMC(Ls{3},options);
    Fd = out1.F(end) + out2.F(end) + out3.F(end);
    varargout = {[varargout,{out1},{out2},{out3}]};
end

p = 1./(1+exp(Fd-Fe));
h = p<.05;

end