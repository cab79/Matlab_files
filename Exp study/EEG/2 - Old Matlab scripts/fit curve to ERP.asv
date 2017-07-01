clear all

subjects = {'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13' 'H14', 'H15', 'H16', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'OA1', 'OA2', 'OA4', 'OA5', 'OA6', 'OA7', 'OA8', 'OA9', 'OA10', 'OA11', 'OA12','OA13'};
  
fnames={
    '_1_avg_ca.mat';
    '_2_avg_ca.mat';    
    '_3_avg_ca.mat';
    '_4_avg_ca.mat';   
    '_5_avg_ca.mat'; 
    '_6_avg_ca.mat'; 
  };

t=1:688;

for n = 34
    subject = subjects(n);
    subject = char(subject);
for x=1:length(fnames);
    fname = fnames(x);
    fname = char(fname);
    sub_fname= [subject fname];
    load(sub_fname); 
    avg = avg(:,1:688);
    
    for e = 1:size(avg,1)
        p2 = polyfit(1:size(avg,2),avg(e,:),2);
        v2=polyval(p2,1:size(avg,2));  
        avg2(e,:)=avg(e,:)-v2; 
        p3 = polyfit(1:size(avg,2),avg(e,:),3);
        v3=polyval(p3,1:size(avg,2));  
        avg3(e,:)=avg(e,:)-v3; 
    end
    plot(t,-avg2(15,:),'b', t,-avg(15,:),'k')
    pause
    name = [subject '___' num2str(x)];
    topoplot(mean(avg([1:2 4:30 33:64],438:500),2), 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-3 3]), colorbar, title(name)
    pause
    topoplot(mean(avg2([1:2 4:30 33:64],438:500),2), 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-3 3]), colorbar, title(name)
    pause
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




subjects = {'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13' 'H14', 'H15', 'H16', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'OA1', 'OA2', 'OA4', 'OA5', 'OA6', 'OA7', 'OA8', 'OA9', 'OA10', 'OA11', 'OA12','OA13'};
  
fnames={
    '_1_avg_ca.mat';
    '_2_avg_ca.mat';    
    '_3_avg_ca.mat';
    '_4_avg_ca.mat';   
    '_5_avg_ca.mat'; 
    '_6_avg_ca.mat'; 
  };

t=1:688;

for n = 1:length(subjects)
    subject = subjects(n);
    subject = char(subject);
    data1 = zeros(64,688);
    data2 = zeros(64,688);
    data3 = zeros(64,688);
for x=1:length(fnames);
    fname = fnames(x);
    fname = char(fname);
    sub_fname= [subject fname];
    load(sub_fname); 
    avg = avg(:,1:688);
    
    for e = 1:size(avg,1)
        p2 = polyfit(1:size(avg,2),avg(e,:),2);
        v2=polyval(p2,1:size(avg,2));  
        avg2(e,:)=avg(e,:)-v2; 
        p3 = polyfit(1:size(avg,2),avg(e,:),3);
        v3=polyval(p3,1:size(avg,2));  
        avg3(e,:)=avg(e,:)-v3; 
    end
    data1 = data1+avg;
    data2 = data2+avg2;
    data3 = data3+avg3;
end
data1 = data1/6;
data2 = data2/6;
data3 = data3/6;
    plot(t,-data3(15,:),'r',t,-data2(15,:),'b', t,-data1(15,:),'k')
    pause
    name = [subject '___' num2str(x)];
    topoplot(mean(data1([1:2 4:30 33:64],438:500),2), 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-3 3]), colorbar, title(name)
    pause
    topoplot(mean(data2([1:2 4:30 33:64],438:500),2), 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-3 3]), colorbar, title(name)
    pause
    topoplot(mean(data3([1:2 4:30 33:64],438:500),2), 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-3 3]), colorbar, title(name)
    pause

end

%%%%%%%%%%%%%%%%%%%%%%%%
%% Apply correction POLY2


clear all

subjects = {'H14'};

fnames={
    '_1_avg_ca.mat';
    '_2_avg_ca.mat';    
    '_3_avg_ca.mat';
    '_4_avg_ca.mat';  
    '_5_avg_ca.mat';
    '_6_avg_ca.mat'; 
  };

for n = 1:length(subjects)
    subject = subjects(n);
    subject = char(subject);
for x=1:length(fnames);
    fname = fnames(x);
    fname = char(fname);
    sub_fname= [subject fname];
    load(sub_fname); 
    
    for e = 1:size(avg,1)
        p2 = polyfit(1:size(avg,2),avg(e,:),2);
        v2=polyval(p2,1:size(avg,2));  
        avg(e,:)=avg(e,:)-v2; 
    end
    eval(['save ' sub_fname ' avg']);
end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%
%% Apply correction POLY3


clear all

subjects = {'F1';'F2'};

fnames={
    '_1_avg_ca.mat';
    '_2_avg_ca.mat';    
    '_3_avg_ca.mat';
    '_4_avg_ca.mat';   
    '_5_avg_ca.mat';
    '_6_avg_ca.mat'; 
  };

for n = 1:length(subjects)
    subject = subjects(n);
    subject = char(subject);
for x=1:length(fnames);
    fname = fnames(x);
    fname = char(fname);
    sub_fname= [subject fname];
    load(sub_fname); 
    
    for e = 1:size(avg,1)
        p3 = polyfit(1:size(avg,2),avg(e,:),3);
        v3=polyval(p3,1:size(avg,2));  
        avg(e,:)=avg(e,:)-v3; 
    end
    eval(['save ' sub_fname ' avg']);
end
    
end