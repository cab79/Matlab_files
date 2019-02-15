%% Analysis: Perceptual model fitting
% This function approaches model inversion from an empirical Bayes
% perspective (based on VBA_MFX function in VBA toolbox), 
% whereby within-subject priors are iteratively refined and
% matched to the inferred parent population distribution.
%  Note: all subjects must use the same model
% See: http://mbb-team.github.io/VBA-toolbox/wiki/VBA-MFX/

function CORE_condor_fit_step2

dbstop if error

pth = pwd;

% add toolbox paths
addpath(genpath(fullfile(pth, 'dependencies')))
addpath(genpath(fullfile(pth, 'Data')))

% get input and output file info
CORE_condor_monitor_outputs;

nsub = importdata('nsub.txt');
nmods = importdata('nmods.txt');

if length(nmods) ~= 2
    error('script assumes there are both perceptual and response models that vary')
end

for pm = 1:nmods(1)
    for rm = 1:nmods(2)
        D_fit=struct;
        
        for d = 1:nsub
            
            % input file index
            ii = (pm-1)*nmods(2)*nsub + (rm-1)*nsub + d;

            % create D_fit
            in_fname = ['input' num2str(ii-1) '.mat'];
            out_fname = ['output' num2str(ii-1) '.mat'];
            load(in_fname,'u','y','S');
            D_fit(d).HGF.u=u;
            D_fit(d).HGF.y=y;
            try
                fit=load(out_fname,'out');
                D_fit(d).HGF.fit=fit.out;
            end
        end

        % save D_fit,S
        save_dir = fullfile(pwd,'Data/fitted');
        if ~exist(save_dir)
            mkdir(save_dir)
        end
        save(fullfile(save_dir,['D_fit_pm' num2str(pm) '_rm' num2str(rm)]),'D_fit','S')
    end
end
quit

function [in,out] = CORE_condor_monitor_outputs

nowtime = now;
in = dir('input*.mat');

fin=0;
while fin==0
    pause(10)
    
    % outputs
    out = dir('output*.mat');
    if length(out)==length(in)
        complete = [out(:).bytes]>0; %& [out(:).datenum]>nowtime;
        if all(complete)
            fin=1;
        elseif sum(complete==0) < 10
            % identify and print stragglers
            filenames = {out(~complete).name};
            disp(['outputs incomplete: ' filenames{:}])
        end
        disp(['number of outputs complete: ' num2str(sum(complete)) '/' num2str(length(complete))])
    end
end

% errors
out = dir('*.err');
if ~isempty(out)
    fsize = [out(:).bytes];
    err_ind = find(fsize>0);
    for n = 1:length(err_ind)
        dat=importdata(out(err_ind(n)).name);
        dat(~cellfun(@isempty,dat))
    end
end

function which_stragglers
% run this code manually
stragglers = [325 524]+1;
nmods = [3 6];
nsub = 44;
strag=struct;
i=0;
for pm = 1:nmods(1)
    for rm = 1:nmods(2)
        for d = 1:nsub
            % input file index
            ii = (pm-1)*nmods(2)*nsub + (rm-1)*nsub + d;
            if any(stragglers==ii)
                i=i+1;
                strag(i).pm = pm;
                strag(i).rm = rm;
                strag(i).d = d;
            end
        end
    end
end