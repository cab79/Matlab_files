clear all
close all

cd('C:\EEGdata');

SUBJECTS = [4:20];
SUBJECTS_n = length(SUBJECTS);


%% RT analysis

ALL_RT=zeros(SUBJECTS_n,11);
ALL_SD=zeros(SUBJECTS_n,11);

for kk=1:SUBJECTS_n
    
    sub=SUBJECTS(1,kk);

% RT=xlsread(['H' num2str(sub) '_pib_orig_rt.xls']); %for healthy subjects
RT=xlsread(['P' num2str(sub) '_pib_orig_rt.xls']); %for patients


RT_n=length(RT);

RT_data=zeros(4,10);
% Columns: fingers from 1 to 10 
% Row 1: count of correct responses
% Row 2: count of wrong responses
% Row 3: percentage of correct responses
% Row 4: percentage of wrong responses


for finger=1:10
    
finger_RT=0;
cor=1;

for i=1:RT_n
    
    if RT(i,1)==finger 
    
        finger_RT(cor)=RT(i,2);
        cor=cor+1;
        
    end 
       
end

finger_RT_mean=mean(finger_RT);
RT_data(1,finger)=finger_RT_mean;

finger_RT_SD=std(finger_RT);
RT_data(2,finger)=finger_RT_SD;

end

ALL_RT(kk,1:10)=RT_data(1,:);
ALL_RT(kk,11)=sub;

ALL_SD(kk,1:10)=RT_data(2,:);
ALL_SD(kk,11)=sub;

end

xlswrite('ALL_Patients_RT.xls',ALL_RT);
xlswrite('ALL_Patients_RT_SD.xls',ALL_SD);
