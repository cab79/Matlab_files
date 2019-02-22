function CORE_eeg_trial_statistics_condor_compile_outputs

pth='/condor_data/cab79/CORE_EEG_stats/';

% add toolbox paths
addpath(genpath(fullfile(pth, 'dependencies_supp')))
addpath(genpath(fullfile(pth, 'Data')))

% get input and output file info
CORE_condor_monitor_outputs;

% get file info
load output0

% S.data_info.D
% S.data_info.C
% S.data_info.CON
% S.data_info.n_chunks

for d = 1:S.data_info.D % subject
    for c = 1:S.data_info.C % channel/component (if used)
        
        if isfield(S,'analysis_type') && strcmp(S.analysis_type,'two_sample') 
        end
            
        % non-parametric Spearman's correlation
        if isfield(S,'analysis_type') && strcmp(S.analysis_type,'SC')
        end
            
        % Multiple regression (linear, non-robust to outliers, non-robust to collinearity)
        if isfield(S,'analysis_type') && strcmp(S.analysis_type,'MR')
        end
            
        % Bayesian regression (linear, non-robust to outliers, non-robust to collinearity)
        if isfield(S,'analysis_type') && strcmp(S.analysis_type,'PEB')
        end
        
        % Ridge regression
        if isfield(S,'analysis_type') && strcmp(S.analysis_type,'RR')
        end
            
        % Bayesian regularised regression
        if isfield(S,'analysis_type') && strcmp(S.analysis_type,'BRR')
            
            for con = 1:S.data_info.CON
                
                for nc = 1:S.data_info.n_chunks
                    
                    % index
                    condor_index = (d-1)*S.data_info.C*S.data_info.CON*S.data_info.n_chunks +(c-1)*S.data_info.CON*S.data_info.n_chunks +(con-1)*S.data_info.n_chunks +nc;
                    index_length = S.data_info.D*S.data_info.C*S.data_info.CON*S.data_info.n_chunks;

                    disp(['compiling output from file ' num2str(condor_index) '/' num2str(index_length)])
                    filein = load(['output' num2str(condor_index-1) '.mat']);
                    out(filein.S.data_info.chunk_index) = filein.out;
                end
                data_dim = filein.S.data_info.dim;

                % reshape
                stats.BRR.alldata(con).b{d,c} = reshape([out(:).muB]',data_dim(1),data_dim(2),[]);
                stats.BRR.alldata(con).s{d,c} = reshape([out(:).muSigma2],data_dim(1),data_dim(2),[]);
                stats.BRR.alldata(con).waic{d,c} = reshape([out(:).waic],data_dim(1),data_dim(2),[]);
                stats.BRR.alldata(con).logl{d,c} = reshape([out(:).logl],data_dim(1),data_dim(2),[]);
                stats.BRR.alldata(con).r2{d,c} = reshape([out(:).r2],data_dim(1),data_dim(2),[]);
                stats.BRR.alldata(con).neglike{d,c} = reshape([out(:).neglike],data_dim(1),data_dim(2),[]);
                stats.BRR.alldata(con).r2test{d,c} = reshape([out(:).r2test],data_dim(1),data_dim(2),[]);
                stats.BRR.alldata(con).xmean{d,c} = reshape([out(:).xmean],data_dim(1),data_dim(2),[]);
                stats.BRR.alldata(con).ymean{d,c} = reshape([out(:).ymean],data_dim(1),data_dim(2),[]);
                stats.BRR.alldata(con).xstd{d,c} = reshape([out(:).xstd],data_dim(1),data_dim(2),[]);
                stats.BRR.alldata(con).ystd{d,c} = reshape([out(:).ystd],data_dim(1),data_dim(2),[]);
                stats.BRR.alldata(con).pred{d,c}=filein.X;
                stats.BRR.alldata(con).pred_group{d,c}=filein.PG;
                stats.trialinfo{con}.idx{d,c}=filein.stats.trialinfo{con}.idx{d,c};
%                 stats.BRR.alldata(con).skew{d,c}=reshape(skew,data_dim(1),data_dim(2),[]);
%                 stats.BRR.alldata(con).kurt{d,c}=reshape(kurt,data_dim(1),data_dim(2),[]);
%                 stats.BRR.alldata(con).hnorm{d,c}=reshape(hnorm,data_dim(1),data_dim(2),[]);
            end   
        end
        
    
    end
end

try
    save(fullfile(S.path.stats,['stats_' S.analysis_type '_' S.data_type '_' S.pred_type{:} '_' S.transform '_' S.sname '.mat']),'stats','S');
    % cleanup
    %delete('input*.mat');
    %delete('output*.mat');
catch
    error('cannot save results')
end
quit


function [in,out] = CORE_condor_monitor_outputs

nowtime = now;
in = dir('input*.mat');

fin=0;
while fin==0
    pause(10)
    out = dir('output*.mat');
    if length(out)==length(in)
        complete = [out(:).bytes]>0; %& [out(:).datenum]>nowtime;
        if all(complete)
            fin=1;
        end
        disp(['number of outputs complete: ' num2str(sum(complete)) '/' num2str(length(complete))])
    end
end