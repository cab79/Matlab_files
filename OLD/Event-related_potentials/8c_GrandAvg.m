clear all

subjects = {'S1_';'S4_';'S5_';'S8_';'S9_';'S10_';'S11_'};
scans = {'pain','nonpain'};

for sess = 1:length(scans)
    
   
    for dt = 1:3
        
        grand_avg=zeros(62,1501);
    
        for sub = 7%1:length(subjects)
            subject = subjects(sub);
            subject = subject{:};
          

            load([subject 'avg_' scans{sess} '_' num2str(dt) '.mat']);
            grand_avg=grand_avg + avg;
            clear avg       

            grand_avg = grand_avg / sub; 

            eval(['save Grand_avg_' scans{sess} '_' num2str(dt) '.mat'  ' grand_avg']);
        end
        clear grand_avg
    end
end




