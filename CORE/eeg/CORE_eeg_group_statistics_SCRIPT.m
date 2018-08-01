clear all
close all
spath = 'C:\Data\CORE\eeg\ana\stats';
sfiles = {
    'stats_SC_all_chan_RT_notrans_20180720T154154.mat' % absolute
    'stats_SC_all_chan_RT_arcsinh_20180720T150942.mat' % absolute
    'stats_SC_all_chan_RT_notrans_20180720T174428.mat' % relative mismatch
    'stats_SC_all_chan_RT_arcsinh_20180720T175126.mat' % relative mismatch
    'stats_SC_all_chan_RT_notrans_20180720T195507.mat' % mismatch, no smoothing/dsample
    'stats_MR_all_chan_RT_notrans_20180720T231023.mat'
    'stats_MR_all_chan_RT_arcsinh_20180720T225208.mat'
    'stats_PEB_all_chan_RT_arcsinh_20180720T222913.mat'
    'stats_RR_all_chan_RT_arcsinh_20180721T081753.mat'
    'stats_BRR_all_chan_RT_arcsinh_20180721T083416.mat' % 100 samples
    'stats_SC_comp_recon_RT_arcsinh_20180721T221137.mat'
    'stats_BRR_all_chan_RT_arcsinh_20180721T172647' % 1000 samples
    };
statfield = 'b';
%statfield = 'rho';
clims = [-5 5]; % t value limits for plotting
xticks = 0:4:600;
param=1;
load('C:\Data\CORE\eeg\ana\prep\chanlocs.mat')

% get data
allstat={};
i=0; % index for fig handles
for f = 12%1:length(sfiles)
    load(fullfile(spath,sfiles{f}));
    allstat{f} = stats;
    
    % groups stats for each analysis
    an = fieldnames(allstat{f});
    for a = 1:length(an)
        da = fieldnames(allstat{f}.(an{a}));
        nc = size(allstat{f}.(an{a}).(da{1}).(statfield),2);
        for c = 1:nc 
            figure('Name',['file_' num2str(f) '_' an{a} ', comp: ' num2str(c)]);pl = 0; % plot index
            for d = 1:length(da)
                disp(['file' num2str(f) ', analysis: ' an{a} ', data type: ' da{d}, ', comp: ' num2str(c)])

                nr = size(allstat{f}.(an{a}).(da{d}).(statfield),1);
            
                clear alldat h p ci t
                for r=1:nr
                    dat=allstat{f}.(an{a}).(da{d}).(statfield){r,c};
                    if ndims(dat)==3
                        dat=dat(:,:,param);
                    elseif ndims(dat)==2 && strcmp(da{d},'GFP')
                        dat=dat(:,param);
                    else
                        dat = dat;
                    end
                    sizdat=size(dat);
                    if ~any(sizdat==1) % if a matrix
                        dat = reshape(dat,prod(sizdat),[]);
                    end
                    alldat(r,:)=dat;
                end
                for s=1:size(alldat,2)
                    [h(s),p(s),~,stt] = ttest(double(alldat(:,s)));
                    t(s)=stt.tstat;
                end
                
                % FDR correction
                [~,fdr_mask] = fdr(p,0.05);
                fdr_p=p.*double(fdr_mask);
                fdr_t=t.*double(fdr_mask);
                % reshape
                if ~any(sizdat==1) 
                    h=reshape(h,sizdat(1),sizdat(2));
                    p=reshape(p,sizdat(1),sizdat(2));
                    t=reshape(t,sizdat(1),sizdat(2));
                    fdr_t=reshape(fdr_t,sizdat(1),sizdat(2));
                    fdr_p=reshape(fdr_p,sizdat(1),sizdat(2));
                end
                % output
                allstat{f}.(an{a}).(da{d}).([statfield '_grpttest']).h=h;
                allstat{f}.(an{a}).(da{d}).([statfield '_grpttest']).p=p;
                allstat{f}.(an{a}).(da{d}).([statfield '_grpttest']).t=t;
                allstat{f}.(an{a}).(da{d}).([statfield '_grpttest']).fdr_t=fdr_t;
                allstat{f}.(an{a}).(da{d}).([statfield '_grpttest']).fdr_p=fdr_p;
                
                hold on
                pl=pl+1;
                subplot(length(da)*2,2,pl)
                imagesc(xticks,[],t,clims); 
                [~,mi]=max(std(t,[],1));
                line(xticks([mi mi]),[1 92],'color','k','linewidth',2)
                title([da{d} ', t values']);
                colorbar
                
                pl=pl+1;
                if ~any(sizdat==1) 
                    subplot(length(da)*2,2,pl)
                    topoplot(t(:,mi),chanlocs);
                end
                
                pl=pl+1;
                subplot(length(da)*2,2,pl)
                imagesc(xticks,[],fdr_t,clims); 
                [~,mi]=max(std(fdr_t,[],1));
                line(xticks([mi mi]),[1 92],'color','k','linewidth',2)
                title([da{d} ', fdr thresholded t values']);
                colorbar
                
                pl=pl+1;
                if ~any(sizdat==1) 
                    subplot(length(da)*2,2,pl)
                    topoplot(fdr_t(:,mi),chanlocs);
                end
            end
        end
    end
    
end

