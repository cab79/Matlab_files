clear all
close all

cd('C:\EEGdata');

SUBJECTS = [4 7:8 10 12:21];
SUBJECTS_n = length(SUBJECTS);


%% Accuracy analysis

ALL_ACC=zeros(SUBJECTS_n,11);

for kk=1:SUBJECTS_n
    
    sub=SUBJECTS(1,kk);

% ACC=xlsread(['H' num2str(sub) '_Exp1_accu.xls']); %for healthy subjects
ACC=xlsread(['P' num2str(sub) '_Exp1_accu.xls']); %for patients


ACC_n=length(ACC);

ACC_data=zeros(4,10);
% Columns: fingers from 1 to 10 
% Row 1: count of correct responses
% Row 2: count of wrong responses
% Row 3: percentage of correct responses
% Row 4: percentage of wrong responses


for finger=1:10
    
finger_correct=0;
cor=1;
finger_wrong=0;
wrong=1;

for i=1:ACC_n
    
    if ACC(i,1)==finger && ACC(i,2)==1
    
        finger_correct(cor)=1;
        cor=cor+1;
        
    elseif ACC(i,1)==finger && ACC(i,2)==0
        
        finger_wrong(wrong)=1;
        wrong=wrong+1;
    
    end 
    
    
end

all_correct=sum(finger_correct);
all_wrong=sum(finger_wrong);
ACC_data(1,finger)=all_correct;
ACC_data(2,finger)=all_wrong;

all=all_correct+all_wrong;
perc_correct=all_correct*100/all;
ACC_data(3,finger)=perc_correct;
ACC_data(4,finger)=100-perc_correct;

end

ALL_ACC(kk,1:10)=ACC_data(3,:);
ALL_ACC(kk,11)=sub;

end

xlswrite('ALL_patients_accuracy.xls',ALL_ACC);
