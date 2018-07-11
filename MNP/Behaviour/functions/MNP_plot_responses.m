function [S,D]= MNP_plot_responses(S,D)
%close all
% PLOT sequence
% if isfield(D(1).Sequence,'adapttype');
%     figure;plot(D(1).Sequence.condnum);
%     hold on; plot(D(1).Sequence.adapttype,'r')
% end
if ismember(S.version,[2 3]);
%     figure;plot(D(1).Sequence.condnum)
%     hold on; 
%     low = ismember(D(1).Sequence.signal(2,:),[1 2]);
%     condplot = nan(1,length(D(1).Sequence.condnum));
%     condplot(low) = D(1).Sequence.condnum(low);
%     s1 = scatter(1:length(D(1).Sequence.condnum),condplot,'k');
%     hold on; 
%     high = ismember(D(1).Sequence.signal(2,:),[3 4]);
%     condplot = nan(1,length(D(1).Sequence.condnum));
%     condplot(high) = D(1).Sequence.condnum(high);
%     s2 = scatter(1:length(D(1).Sequence.condnum),condplot,'r');
%     title('design')
%     legend([s1,s2],'low discrim','high discrim')
%     
    figure;plot(D(1).Sequence.condnum)
    hold on; 
    low = ismember(D(1).Sequence.signal(2,:),[1 3]);
    condplot = nan(1,length(D(1).Sequence.condnum));
    condplot(low) = D(1).Sequence.condnum(low);
    s1 = scatter(1:length(D(1).Sequence.condnum),condplot,'k');
    hold on; 
    high = ismember(D(1).Sequence.signal(2,:),[2 4]);
    condplot = nan(1,length(D(1).Sequence.condnum));
    condplot(high) = D(1).Sequence.condnum(high);
    s2 = scatter(1:length(D(1).Sequence.condnum),condplot,'r');
    title('design')
    legend([s1,s2],'low intensity','high intensity')
    
    figure;plot(D(1).Sequence.condnum)
    hold on; 
    low = ismember(D(1).Sequence.signal(1,:),1);
    condplot = nan(1,length(D(1).Sequence.condnum));
    condplot(low) = D(1).Sequence.condnum(low);
    s1 = scatter(1:length(D(1).Sequence.condnum),condplot,'k');
    hold on; 
    high = ismember(D(1).Sequence.signal(1,:),2);
    condplot = nan(1,length(D(1).Sequence.condnum));
    condplot(high) = D(1).Sequence.condnum(high);
    s2 = scatter(1:length(D(1).Sequence.condnum),condplot,'r');
    title('design')
    legend([s1,s2],'cue type 1','cue type 2')
end

%
if S.accuracy.on

    %5. PLOT %correct as a moving average over time
    for d = 1:length(D)
        Y=D(d).Processed.macorrect;
        figure
        plot(Y,'linewidth',2)
        ylabel('% correct')
        xlabel('trials')
        title(['% correct, moving average of ' num2str(S.movingavg) ' trials'])
    end

    %5. PLOT difference in %correct over time between valid and invalid cues for each block 
    for d = 1:length(D)
        figure
        ii=0;
        for b = 1:length(D.Processed.blockcorrectmovavg{1})
            plotind = find(~isnan(D.Processed.blockcorrectmovavg{1}{b}));
            if ~isempty(plotind) && any(D.Processed.blockcorrectmovavg{1}{b}(plotind)) && length(plotind)>S.movingavg
                ii=ii+1;
                subplot(2,2,ii)
                Y=D.Processed.blockcorrectmovavg{1}{b};
                Yind = 1:length(Y);
                % clear leading zeros
                index = find(Y ~= 0, 1, 'first');
                Y = Y(index:end);
                Yind = Yind(index:end);
                % find values less than an equal/greater than S.movingavg and
                % plot separately
                hold on
                Yplot = Y; Yplot(Yind<S.movingavg)=NaN;
                scatter(1:length(Yplot),Yplot,'b','filled')
                lsline
                Yplot = Y; Yplot(Yind>=S.movingavg)=NaN;
                scatter(1:length(Yplot),Yplot,'b')
                hold off
                ylabel('% correct, valid-invalid')
                xlabel('trials')
                title(['Block ' num2str(b)])
            end
        end
    end

    %5. PLOT %correct
    for d = 1:length(D)
        for tm = 1:length(S.trialmax)
            Y=D(d).Processed.condcorrectfract{tm};
            X=1:length(Y);
            figure
            bar(X,Y)
            labels = arrayfun(@(value) num2str(value,'%2.0f'),cell2mat(D(d).Processed.numtrials{tm}),'UniformOutput',false);
            text(X,Y,labels,'HorizontalAlignment','center','VerticalAlignment','bottom') 
              % clears X axis data
              %set(gca,'XTick',[]);
            ylabel('fraction correct')
            xlabel('condition')
            switch S.version
                case 1
                    %labels = {'adaptive: low','adaptive: high','low prob: low','low prob: high','equal prob: low','equal prob: high','high prob: low','high prob: high',};
                    labels = {'low prob: low','low prob: high','equal prob: low','equal prob: high','high prob: low','high prob: high',}; % v1
                case {2 3}
                    labels = {'high prob: pair 1','low prob: pair 2','equal prob: pair 1','equal prob: pair 2','high prob: pair 2','low prob: pair 1',}; % v2
            end
            set(gca,'xticklabel', labels)
            %title([D(d).subname ', trial number: ' num2str(S.trialmax{tm})])
            title([D(d).subname])
            hold on
            plot(xlim,[0.5 0.5], 'k--')
        end
    end

    %5. PLOT %correct for some conditions only, across blocks
    %close all
    if S.version==1
        plotcond = [3 4];
        for d = 1:length(D)
            for tm = 1:length(S.trialmax)
                figure
                dat=cell2mat(D(d).Processed.blockcondcorrectfract{tm}(plotcond)')';
                b = bar([dat])
                ylabel('fraction correct')
                xlabel('block')
                legend(b,{'low','high'});
                title(D(d).subname)
                hold on
                plot(xlim,[0.5 0.5], 'k--')
            end
        end
    end

    % PLOT by stimulus
    if isfield(D.Processed,'stimcondcorrectfract')
        for d = 1:length(D)
            for tm = 1:length(S.trialmax)

                % for each condition, breaking down by stimulus type
                Y=cell2mat(D.Processed.stimcondcorrectfract{:}');
                figure
                b=bar(Y)
                ylabel('fraction correct')
                xlabel('condition')
                switch S.version
                    case 1
                        %labels = {'adaptive: low','adaptive: high','low prob: low','low prob: high','equal prob: low','equal prob: high','high prob: low','high prob: high',};
                        labels = {'low prob: low','low prob: high','equal prob: low','equal prob: high','high prob: low','high prob: high',}; % v1
                    case {2 3}
                        labels = {'high prob: pair 1','low prob: pair 2','equal prob: pair 1','equal prob: pair 2','high prob: pair 2','low prob: pair 1',}; % v2
                end
                set(gca,'xticklabel', labels)
                %title([D(d).subname ', trial number: ' num2str(S.trialmax{tm})])
                title([D(d).subname])
                legend(b,{'low','high','low x2','high x2'});
                hold on
                plot(xlim,[0.5 0.5], 'k--')

                % plot stims averaged across all conditions
                Y = nanmean(Y,1);
                figure
                b=bar(Y)
                ylabel('fraction correct')
                xlabel('stimulus type')
                switch S.version
                    case {2 3}
                        labels = {'low','high','low x2','high x2'}; % v2
                end
                set(gca,'xticklabel', labels)
                %title([D(d).subname ', trial number: ' num2str(S.trialmax{tm})])
                title([D(d).subname])
                hold on
                plot(xlim,[0.5 0.5], 'k--')
            end
        end
    end

    % PLOT by cue/stimulus
    if isfield(D.Processed,'stimcuecorrectfract')
        for d = 1:length(D)
            for tm = 1:length(S.trialmax)

                % for each cue, breaking down by stimulus type
                Y=cell2mat(D.Processed.stimcuecorrectfract{:}');
                figure
                b=bar(Y)
                ylabel('fraction correct')
                xlabel('cue type')
                switch S.version
                    case {2 3}
                        labels = {'high pitch','low pitch',}; % v2
                end
                set(gca,'xticklabel', labels)
                %title([D(d).subname ', trial number: ' num2str(S.trialmax{tm})])
                title([D(d).subname])
                legend(b,{'low','high','low x2','high x2'});
                hold on
                plot(xlim,[0.5 0.5], 'k--')

            end
        end
    end
end

if S.RT.on
    %5. PLOT %correct as a moving average over time
    for d = 1:length(D)
        Y=D(d).Processed.malogrt;
        figure
        plot(Y,'linewidth',2)
        ylabel('log(RT)')
        xlabel('trials')
        title(['log(RT), moving average of ' num2str(S.movingavg) ' trials'])
    end
    
    %5. PLOT %correct
    for d = 1:length(D)
        for tm = 1:length(S.trialmax)
            Y=D(d).Processed.cond_logrt{tm};
            X=1:length(Y);
            figure
            bar(X,Y)
            %labels = arrayfun(@(value) num2str(value,'%2.0f'),cell2mat(D(d).Processed.numtrials{tm}),'UniformOutput',false);
            %text(X,Y,labels,'HorizontalAlignment','center','VerticalAlignment','bottom') 
            ylabel('log(RT)')
            xlabel('condition')
            switch S.version
                case 1
                    %labels = {'adaptive: low','adaptive: high','low prob: low','low prob: high','equal prob: low','equal prob: high','high prob: low','high prob: high',};
                    labels = {'low prob: low','low prob: high','equal prob: low','equal prob: high','high prob: low','high prob: high',}; % v1
                case {2 3}
                    labels = {'high prob: pair 1','low prob: pair 2','equal prob: pair 1','equal prob: pair 2','high prob: pair 2','low prob: pair 1',}; % v2
            end
            set(gca,'xticklabel', labels)
            %title([D(d).subname ', trial number: ' num2str(S.trialmax{tm})])
            title([D(d).subname])
        end
    end
    
    % PLOT by stimulus
    if isfield(D.Processed,'stimcond_logrt')
        for d = 1:length(D)
            for tm = 1:length(S.trialmax)

                % for each condition, breaking down by stimulus type
                Y=cell2mat(D.Processed.stimcond_logrt{:}');
                figure
                b=bar(Y)
                ylabel('log(RT)')
                xlabel('condition')
                switch S.version
                    case 1
                        %labels = {'adaptive: low','adaptive: high','low prob: low','low prob: high','equal prob: low','equal prob: high','high prob: low','high prob: high',};
                        labels = {'low prob: low','low prob: high','equal prob: low','equal prob: high','high prob: low','high prob: high',}; % v1
                    case {2 3}
                        labels = {'high prob: pair 1','low prob: pair 2','equal prob: pair 1','equal prob: pair 2','high prob: pair 2','low prob: pair 1',}; % v2
                end
                set(gca,'xticklabel', labels)
                %title([D(d).subname ', trial number: ' num2str(S.trialmax{tm})])
                title([D(d).subname])
                legend(b,{'low','high','low x2','high x2'});

                % plot stims averaged across all conditions
                Y = nanmean(Y,1);
                figure
                b=bar(Y)
                ylabel('log(RT)')
                xlabel('stimulus type')
                switch S.version
                    case {2 3}
                        labels = {'low','high','low x2','high x2'}; % v2
                end
                set(gca,'xticklabel', labels)
                %title([D(d).subname ', trial number: ' num2str(S.trialmax{tm})])
                title([D(d).subname])
            end
        end
    end
end

% PLOT adaptive thresholds
if isfield(D(1).Sequence,'adapttype');
    n_rev = 6;
    av_type = 1;
    av_para = []; col = {'r','m','y'};
    for d = 1:length(D)
        for atype = 1:2
            ind = find(D(d).Output.adaptive(:,10)==atype & ~isnan(D(d).Output.adaptive(:,7)));
            if ~isempty(ind)
                thresh = D(d).Output.adaptive(ind,7);
                rev = D(d).Output.adaptive(ind,3);
                block = D(d).Output.adaptive(ind,7);

                % re-calc thresholds
                %if n_rev
                %    thresh = nan(length(thresh),1);
                %    for i = n_rev:length(thresh)
                %        thresh(i)=mean(rev(i-(n_rev-1):i,1));
                %    end
                %end

                av_thresh = nan(length(av_para),length(thresh));
                switch av_type
                    case 1
                        % moving average
                        for av = 1:length(av_para)
                            for i = av_para(av):length(thresh)
                                av_thresh(av,i) = mean(thresh((i-av_para(av)+1):i,1));
                            end
                        end
                    case 2
                        % polynomial
                        for av = 1:length(av_para)
                            p = polyfit(1:length(thresh),thresh',av_para(av));
                            av_thresh(av,:)=1:length(thresh); 
                            av_thresh(av,:)=av_thresh(av,:)*p(1)+p(2);
                            av_thresh(av,:)=av_thresh(av,:)';
                        end
                    case 3
                        % smooth
                        for av = 1:length(av_para)
                            for i = 1:length(thresh)
                                sm_thresh = smooth(thresh(1:i,1),min(i,av_para(av)),'lowess');
                                av_thresh(av,i) = sm_thresh(end);
                            end
                        end
                    case 4
                        % weighted moving average 
                        for av = 1:length(av_para)
                            for i = av_para(av):length(thresh)
                                %alpha = 0.06;
                                %sm_thresh = filter(alpha, [1 alpha-1], thresh((i-av_para+1):i,1));
                                W = 1:av_para(av); % relative weights
                                sm_thresh = conv(thresh((i-av_para(av)+1):i,1)', W./sum(W), 'full'); % linear weighting
                                av_thresh(av,i) = sm_thresh(av_para(av));
                            end
                        end
                end

                % identify stabilty
                stab_type = 3;
                maxpercentchange = 0.5;
                npairs = 20;


                stable_thresh = nan(av,length(thresh));
                if ismember(stab_type,[1 2]);
                    for av = 1:length(av_para)
                        mavg = av_thresh(av,:);
                        mavg_ind = find(~isnan(mavg));
                        range = max(thresh)-min(thresh);
                        for m = npairs+1:length(mavg_ind)
                            switch stab_type
                                case 1
                                    % linear fit
                                    p = polyfit(1:npairs,mavg(mavg_ind(m-npairs+1):mavg_ind(m)),1);
                                    slope=abs(p(1)); 
                                    if (slope/range)*100<maxpercentchange
                                        stable_thresh(av,mavg_ind(m)) = mavg(mavg_ind(m));
                                    end

                                case 2
                                    % 
                                    diffval = [];
                                    for i = 1:npairs
                                        diffval(i) = abs(diff([mavg(mavg_ind(m-(i-1)-1)),mavg(mavg_ind(m-(i-1)))]));
                                    end
                                    if all((diffval/range)*100<maxpercentchange)
                                        stable_thresh(av,mavg_ind(m)) = mavg(mavg_ind(m));
                                    end

                            end
                        end
                    end
                elseif ismember(stab_type,3);
                    % get max mavg length and it's index
                    [maxlen,maxi] = max(av_para);
                    % indices of max mavg
                    maxmavg = av_thresh(maxi,:);
                    mavg_ind = find(~isnan(maxmavg));
                    for m = 1:length(mavg_ind)
                        trend = [];
                        for av = 1:length(av_para)-1
                            trend(av) = (av_thresh(av,mavg_ind(m))-av_thresh(av+1,mavg_ind(m)));
                        end
                        if ~all(trend>0) && ~all(trend<0)
                            stable_thresh(1,mavg_ind(m)) = mean(av_thresh(:,mavg_ind(m)));
                        end
                    end
                end
                figure
                scatter(1:length(thresh),thresh,'b');
                hold on
                for av = 1:length(av_para)
                    scatter(1:length(thresh),av_thresh(av,:),col{av});
                end
                %for av = 1:length(av_para)
                %    scatter(1:length(thresh),stable_thresh(av,:),'k','filled');
                %end
                scatter(1:length(thresh),D(d).Output.adaptive(ind,11),'k','filled');
                highval=D(d).Output.adaptive(end,11)+2*std(thresh);
                title([num2str(D(d).Output.adaptive(end,11)) ', high: ' num2str(highval)])
                hold off
                figure
                scatter(1:length(thresh)-1,diff(thresh),'b');
            end

        end
    end
end

% plot RT distribution
for d = 1:length(D)
    figure
    hist(D(d).Output.RT);
    title(D(d).subname)
end