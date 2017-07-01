clear all
close all

filepath = 'C:\Data\CRPS-DP\Behaviour - accuracy data';
run('M:\Matlab\Matlab_files\CRPS_digits\loadsubj.m');

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
        
        rep_idx = find(ACC(:,1)==finger & ~isnan(ACC(:,2)));
        ACC(rep_idx,2) = res;
        for i = rep_idx
            ACC(i,3)=ACC(i,1)==ACC(i,2);
        end
        
      
    end
    
    [pth nme ext] = fileparts(fname);
    fname_save = [nme '_corrected.csv'];
    csvwrite(fname_save,ACC);
   
end


