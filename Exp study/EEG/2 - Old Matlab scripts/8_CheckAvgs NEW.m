

% Topoplots: grand averages
clear all
chan_locs = 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan_MR62.locs';
ele = [5,6,17,18,19]; % electrodes to average when searching for peaks
for dt = 1:2
    load(['Grand_avg_' num2str(dt) '_ca.mat'])
    subplot(2,4,(dt-1)*4+1); 
    [handle,Zi,grid] = topoplot(mean(grand_avg(:,750:1000),2), chan_locs, 'maplimits','absmax'); title(['early_' num2str(dt)]);
    hold on
    subplot(2,4,(dt-1)*4+2); 
    [handle,Zi,grid] = topoplot(mean(grand_avg(:,1750:2000),2), chan_locs, 'maplimits','absmax'); title(['late_' num2str(dt)]);
    hold on
    subplot(2,4,(dt-1)*4+3); 
    [XMAX,IMAX,XMIN,IMIN] = extrema(-mean(grand_avg(ele,2050:2300),1));%%search for extrema on internet - peak detection prog, xmax ampli imax latency
    LATn = IMAX(1)+2049;
    [handle,Zi,grid] = topoplot(mean(grand_avg(:,LATn-10:LATn+10),2), chan_locs, 'maplimits','absmax'); title(['n2_' num2str(dt)]);
    hold on
    subplot(2,4,(dt-1)*4+4); 
    [XMAX,IMAX,XMIN,IMIN] = extrema(mean(grand_avg(ele,2100:2500),1));%%search for extrema on internet - peak detection prog, xmax ampli imax latency
    LAT = IMAX(1)+2099;
    [handle,Zi,grid] = topoplot(mean(grand_avg(:,LAT-10:LAT+10),2), chan_locs, 'maplimits','absmax'); title(['p2_' num2str(dt)]);
    hold on
end

% Topoplots: individual subjects
clear all
chan_locs = 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan_MR62.locs';
subjects = {'P1_';'P2_';'P5_';'P6_';'P7_';'P8_';'P15_';'P16_';'S1_';'S2_';'S3_';'S5_';'S6_';'S8_';'S9_';'S10_';'S11_';'S18_';'S20_';'S4_';'P14_';'P20_';'P24_';'P30_';'P32_';'P33_';'P19_'};
load N2latencies
load P2latencies
load N1latencies
N2 = (N2LAT+4000)/2;
P2 = (P2LAT+4000)/2;
N1 = (N1LAT+4000)/2;
for sub = 1:length(subjects)
    subject = subjects(sub);
    subject = subject{:};
    for dt = 1:2
        load([subject 'avg_' num2str(dt) '_ca.mat'])
        avg = eegfilt(avg,500,2,20); %filter the data
        avg = blcorrect4(avg, 1750);
        subplot(1,2,dt); 
        lat = N1(sub,dt);
        [handle,Zi,grid] = topoplot(mean(avg(:,lat-5:lat+5),2), chan_locs, 'maplimits','absmax'); title([subject num2str(dt)]);
        hold on;   
    end
    pause
end

% Topoplot loop: grand averages
clear all
chan_locs = 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan_MR62.locs';
load(['Grand_avg_1_ca.mat'])
grand_avg = eegfilt(grand_avg,500,2,20); %filter the data
grand_avg = blcorrect4(grand_avg, 1750);
step = 5;
for i = 2000:step:size(grand_avg,2)
    [handle,Zi,grid] = topoplot(mean(grand_avg(:,i:i+step),2), chan_locs, 'maplimits','absmax'); 
    time = i*2-4000;
    title([num2str(time) 'ms post stimulus'])
    pause
end

% Check Averages
clear all
t=0:2:1498;
col = {'b','r'};
for dt = 1:2
    load(['Grand_avg_' num2str(dt) '_ca.mat'])
    plot(t, -grand_avg(18,t/2+2000), col{dt})
    hold on
end

% Topoplot loop: individual subjects
clear all
subjects = {'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13' 'H14', 'H15', 'H16', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'OA1', 'OA2', 'OA4', 'OA5', 'OA6', 'OA7', 'OA8', 'OA9', 'OA10', 'OA11', 'OA12','OA13','OA14', 'OA15', 'OA16','OA17'};
chan_locs = 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan_MR62.locs';
load N2latencies
load P2latencies
load N1latencies
N2 = (N2LAT+4000)/2;
P2 = (P2LAT+4000)/2;
N1 = (N1LAT+4000)/2;
step = 5;
sub = 1;
subject = subjects(sub);
subject = subject{:};
dat1 = load([subject 'avg_1_ca.mat']);
dat2 = load([subject 'avg_2_ca.mat']);
grand_avg = (dat1.avg+dat2.avg)/2;
grand_avg = eegfilt(grand_avg,500,2,20); %filter the data
grand_avg = blcorrect4(grand_avg, 1750);
for dt = 1:2
    figure(1);
    subplot(1,2,dt); 
    lat = N1(sub,dt);
    eval(['avg = dat' num2str(dt) '.avg;']);
    avg = eegfilt(avg,500,2,20); %filter the data
    avg = blcorrect4(avg, 1750);
    [handle,Zi,grid] = topoplot(mean(avg(:,lat-5:lat+5),2), chan_locs, 'maplimits','absmax'); 
    title([subject num2str(dt) ' ' num2str(lat*2-4000) 'ms post stimulus']);
    hold on;   
end
for i = 2075:step:size(grand_avg,2)
    figure(2);
    [handle,Zi,grid] = topoplot(mean(grand_avg(:,i:i+step),2), chan_locs, 'maplimits','absmax'); 
    time = i*2-4000;
    title([num2str(time) 'ms post stimulus'])
    pause
end


% Check Averages on individual subjects
clear all
close all
t=0:2:1498;
col = {'b','r'};
col2 = {'b--','r--'};
subjects = {'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13' 'H14', 'H15', 'H16', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'OA1', 'OA2', 'OA4', 'OA5', 'OA6', 'OA7', 'OA8', 'OA9', 'OA10', 'OA11', 'OA12','OA13','OA14', 'OA15', 'OA16','OA17'};
for sub = 1:length(subjects)
    figure('units','normalized','outerposition',[0 0 1 1]);
    hold on
    subject = subjects(sub);
    subject = subject{:};
    dti = 0;
    for dt = [4 6]
        dti=dti+1;
        load([subject '_' num2str(dt) '_avg_ca.mat'])
        avg = eegfilt(avg,500,2,20); %filter the data
        avg = blcorrect4(avg, 1750);
        figure(1)
        plot(t, -mean(avg([14 15 20],t/2+2000),1), col{dti}); title(subject);
        hold on
        plot(t, -mean(avg([18 24 25 47 53 57],t/2+2000),1), col2{dti}); title(subject);
        hold on
    
    %check latency
    load N2latencies
    load P2latencies
    load N1latencies
    load amplitudes
    AMP = reshape(AMP2,length(subjects),6,5);
    N2 = (N2LAT+4000)/2;
    P2 = (P2LAT+4000)/2;
    N1 = (N1LAT+4000)/2;
    Y = zeros(1,length(t));
    Z = zeros(1,length(t));
    Q = zeros(1,length(t));

    Y(P2(sub,dt)-2000-t(1)/2) = AMP(sub,dt,5);
    plot (t, -Y(:,:), 'y.')
    hold on


    Z(N2(sub,dt)-2000-t(1)/2) = AMP(sub,dt,4);
    plot (t, -Z(:,:), 'm.')
    hold on

    
    Q(N1(sub,dt)-2000-t(1)/2) = AMP(sub,dt,3);
    plot (t, -Q(:,:), 'g.')
    hold on

    
    end

    pause
    close all
end
