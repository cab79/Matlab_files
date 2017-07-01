%%% boucle electrode par electrode, sample par sample, test-t 
close all;
clear all;



Good_subj=[1];

% 2 7 12 14 19 25 27 28 29

nsubj =length(Good_subj);


for k= 1:nsubj
    
    i=Good_subj(1,k);
    
    
    
    file1=['P20_Exp2.Left.aff.set']; %red(magenta)
    file2=['P20_Exp2.Right.Unaff.set']; %blue(cyan)
    
%     file1=['nogo__' num2str(i) '_filtd_avgd_epochd_rmvbase_itpld_goodchan_beh_epdeleted_ICA_pruned_600_st2winB100.set']; %red
%     file2=['nogo__' num2str(i) '_filtd_avgd_epochd_rmvbase_itpld_goodchan_beh_epdeleted_ICA_pruned_200_st2winB100.set']; %blue
    
    directory='C:\EEGdata/';
    n_elec=92;
    threshold = 0.05;
    temporal_window_threshold = 5; %physical time or sampling point
    
    
    
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    
    fig = figure( 'Name', 'P4: gr10 vs. gr12')
    EEG = pop_loadset(file1, directory);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    EEG = pop_loadset(file2, directory);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    [a,b,c] = size(EEG.data);
    %%% Dans EEG.data, a=nombre d'�lectrodes, b=nombre de points par epoch,
    %%% c=nombre d'epochs
    ttest_elec = zeros(n_elec,b);
    % Statistical thresholds for t-tests
    thresholded_data = zeros(n_elec,b);
    double_thresholded_data =zeros(n_elec,b);
    
    for elec=1:n_elec
        
        [a,b,c] = size(ALLEEG(1).data);
        
        cond1 = zeros(elec,b);
        cond1_trials = zeros(elec,c,b);
        for time = 1:b;
            cond1(elec,time)= mean(ALLEEG(1).data(elec,time,:));
            for trials = 1:c;
                cond1_trials(elec,trials,time)= ALLEEG(1).data(elec,time,trials);
            end
        end
        
        [a,b,c] = size(ALLEEG(2).data)
        
        cond2 = zeros(elec,b);
        cond2_trials = zeros(elec,c,b);
        
        for time = 1:b;
            cond2(elec,time)= mean(ALLEEG(2).data(elec,time,:));
            for trials = 1:c;
                cond2_trials(elec,trials,time)= ALLEEG(2).data(elec,time,trials);
            end
            [t,p] = ttest2(cond2_trials(elec,:,time),cond1_trials(elec,:,time));
            ttest_elec(elec,time) = p;
            thresholded_data(elec,time) = 9999;
            if p <= threshold
                thresholded_data(elec,time) = 0;
            end
        end
        
        %% temporal window threshold :
        for zz = 1:(b-1);
            window = 0;
            if (thresholded_data(elec,zz) == 0)
                for prec = zz:-1:1
                    if (thresholded_data(elec,prec) == 0)
                        window = window + 1;
                    end
                    if (thresholded_data(elec,prec) == 9999)
                        break
                    end
                end
                for subs = (zz+1):b
                    if (thresholded_data(elec,subs) == 0)
                        window = window + 1;
                    end
                    if (thresholded_data(elec,subs) == 9999)
                        break
                    end
                end
            end
            double_thresholded_data(elec,zz) = 9999;
            if (window >= temporal_window_threshold)
                double_thresholded_data(elec,zz) = 0;
            end
        end
        
        
        for aa = 1:250
            %     Before it was: aa = 1:850; %set frames per epoch
            newtime(aa) = aa*4-200; % aa*(1000/sample rate) - baseline - Tristan did it for -200:3000 epochs
            %before it was: newtime(aa) = aa*2-200; % aa*(1000/sample rate) - baseline - Tristan did it for -200:3000 epochs
        end
        
        
        
        lin=fix(sqrt(n_elec));
        subplot (lin+1,lin+1,elec)
        cond3=cond2-cond1;
        %subplot(4,4,elec)
        name=strcat( 'electrode', int2str(elec) ) ;
        plot(newtime(1,1:250),cond2(elec,1:250), 'b','LineWidth',2) % instead of blue
        hold on
        plot(newtime(1,1:250),cond1(elec,1:250),'r','LineWidth',2) % instead of red
        % plot(newtime(1,1:850),cond3(elec,1:850),'k')
        %plot(newtime(1,1:850),(ttest_elec(elec,1:850)),'b')
        plot(newtime(1,1:250),double_thresholded_data(elec,1:250),'.g','LineWidth',2)
        axis([-200 796 -6 6])
        title(name,'fontsize',9)
        h=gca;
        set(h, 'FontSize', 9)
        
    end
    
    
    com = [ 'axis on;' ...
        'clear xlabel ylabel;' 'xlabel(''''Time (ms)'''');' ...
        'ylabel(''''Voltage (\muV)'''');' ];
    
    axcopy(gcf,com)
    
fig;
%     saveas(gcf,['t-test_nogo_' num2str(i) '_600_200_st2win_base100.fig'])
    saveas(gcf,['t_test_P4_gr10_gr12.fig'])
    
end
