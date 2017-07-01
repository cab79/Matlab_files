clear all

subjects = {'S1_';'S4_';'S5_';'S8_';'S9_';'S10_';'S11_';'Average'};
scans = {'pain','nonpain'};

load gamma1.mat 
load gamma2.mat 
load beta.mat 
load alpha.mat 
load theta.mat 
load delta.mat 
chan_locs = 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan_MR62.locs';

data = alpha;
period = [1:3];

data(:,length(subjects),:,:) = mean(data(:,:,:,:),2);

for sub = 1:length(subjects)
    subject = subjects(sub);
    subject = subject{:};

for per = period
    
    subt = zeros(62,1);
    
for sess = 1:length(scans)
    
    %text = (['Subject ' subject(1:length(subject)-1) ', Session ' session ', Period ' num2str(per)]);
    dat = data(:,sub,sess,per);
    subplot(length(period),length(scans)+1,(3*per-2)+(sess-1)); 
    [handle,Zi,grid] = topoplot(dat, chan_locs, 'maplimits','absmax'); 
    if per == 1 && sess ==1
    Zimin = 0.9*min(min(Zi));
    Zimax = 0.9*max(max(Zi));
    end
    caxis([Zimin Zimax])
    colorbar
    
    subt = data(:,sub,sess,per) - subt;
    
    %if data == delta
    %    caxis([-20 20]),
    %elseif data == alpha 
    %    caxis([-20 20]),
    %elseif data == theta
    %    caxis([-20 20]),
    %elseif data == beta
    %    caxis([-10 10]),
    %elseif data == gamma1
    %    caxis([-20 0]),
    %elseif data == gamma2
    %    caxis([-20 0]),
    %end
end
    subplot(length(period),length(scans)+1,(3*per)); 
    [handle,Zi,grid] = topoplot(subt, chan_locs, 'maplimits','absmax'); 
    mtit(['Subject ' subject(1:length(subject)-1)]);

%% Optional inputs:
   % 'maplimits'       - 'absmax'   -> scale map colors to +/- the absolute-max (makes green 0); 
   %                     'maxmin'   -> scale colors to the data range (makes green mid-range); 
   %                     [lo.hi]    -> use user-definined lo/hi limits {default: 'absmax'}
   % 'style'           - 'map'      -> plot colored map only
   %                     'contour'  -> plot contour lines only
   %                     'both'     -> plot both colored map and contour lines
   %                     'fill'     -> plot constant color between contour lines
   %                     'blank'    -> plot electrode locations only {default: 'both'}
   % 'electrodes'      - 'on','off','labels','numbers','ptslabels','ptsnumbers' See Plot detail 
   %                     options below. {default: 'on' -> mark electrode locations with points
   %                     unless more than 64 channels, then 'off'}. 
    
end
pause
end