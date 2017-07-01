clear all
close all

filepath = 'C:\Data\CRPS-DP\Behaviour - accuracy data';
run('C:\Data\Matlab\Matlab_files\CRPS_digits\loadsubj.m');


cd(filepath)
files = dir('*accu.csv');

SUBJECTS_n = length(files);

%% Accuracy analysis

ALL_ACC=cell(SUBJECTS_n,21);

digit_ranges = {sort(1:5,'descend');6:10};

for kk=1:SUBJECTS_n
    
    fname = files(kk).name;
    C = strsplit(fname,'_');
    sub=C{1};

    ACC=xlsread(fname); %for patients
    ACC_n=length(ACC);
    ACC_data=zeros(21,10);
    % Rows 1 to 10: Count of responses from 1 to 10
    % Rows 11 to 20: Percentage of responses from 1 to 10

    for finger=1:10
        res_count = zeros(10,1);
        res = ACC(find(ACC(:,1)==finger),2);
        res(isnan(res)) = [];
        
        for s = 1:size(digit_ranges)
            if any(digit_ranges{s}==finger)
                r = digit_ranges{s}; % find out which fingers are on the hand we are dealing with
                nri = 1:2;
                nri(nri==s)=[];
                nr = digit_ranges{nri};
            end
        end
        for e = 1:length(res)
            if ~any(r==res(e)) && res(e)~=0
                res(e)=r(nr==res(e));
            end
        end
        
        for d=1:10
            res_count(d,1) = sum(res(:) == d);
        end
        res_perc = res_count*100 ./ sum(res_count);
        ACC_data(1:10,finger) = res_count;
        ACC_data(11:20,finger) = res_perc;
        ACC_data(21,finger) = res_perc(finger);
    end
    
    for t=2:length(ACC)
        ACC(t,4) = abs(ACC(t,1)-ACC(t-1,1));
    end
    
    for change=1:5
       chgacc(change,1) = sum(ACC(:,4)==change-1 & ACC(:,3)==1);
       chgall(change,1) = sum(ACC(:,4)==change-1 & ~isnan(ACC(:,2)));
    end
    chgperc = chgacc*100 ./ chgall;
    
    ALL_ACC{kk,1}=sub;
    ALL_ACC(kk,2:11)=num2cell(ACC_data(21,:));
    ALL_ACC(kk,12:111)=num2cell(reshape(ACC_data(11:20,:),1,100));
    ALL_ACC(kk,112:116)=num2cell(chgperc');

end

xlswrite('ALL_patients_accuracy.xls',ALL_ACC);
