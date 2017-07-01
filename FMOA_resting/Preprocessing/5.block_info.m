subjects = {'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13' 'H14', 'H15', 'H16', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'OA1', 'OA2', 'OA4', 'OA5', 'OA6', 'OA7', 'OA8', 'OA9', 'OA10', 'OA11', 'OA12', 'OA13', 'OA14', 'OA15', 'OA16', 'OA17'};

b_order = [1 2 1 2 1 2 1 2 1 2 1 2 1 2 2 1 2 1 2 1 1 2 2 1 2 1 2 1 2 1 2 1 1 2 2 2 1 1 2 1 2 1 2 2 1 2 2 1]; % 1 = open first 

SD_thresh = 10;
SD_min = 2;

exclude_none = 0;
auto_exclude = 0;
manual_exclude = 1;

excludes = cell(length(subjects),1);

if manual_exclude ==1;
    excludes = {
    '[10,53;]'%1
    '[1,20,39;]'%2
    '[32,58;]'%3
    '[16,18,25,26,27,28,29,30,42;]'%4
    '[16,17,18,19,20;]'%5
    ''%6
    '[3,8,9,17;]'%7
    '[2,7,16,18,32,33;]'%8
    '[30,48;]'%9
    '[30,31,49;]'%10
    '[1,16,30,31,44,49;]'%11
    '[30,44,46,49,56,57;]'%12
    '[14,15,30,32,35;]'%13
    '[3,4,5,6,7,16,27,31,50;]'%14
    ''%15
    '[15,30;]'%16
    '[34;]'%1
    '[6,15,30,31,46;]'%2
    '[6,24,31;]'%3
    '[23,33,37,38,43;]'%4
    '[16;]'%5
    '[6,9,10,11,16;]'%6
    '[4,8,16,23,25,42,43,56,57,58;]'%7
    '[1;]'%8
    '[1,2,3;]'%9
    '[7;]'%10
    '[18,20,46;]'%11
    '[15,16;]'%12
    '[27;]'%13
    '[19;]'%14
    '[1,20,23,24,27,28,40,41,42;]'%15
    ''%16
    '[16,17,18,19,20,21,22,23,24,25,26,27,28,29,35,36,37,38,39,40:58;]'%1
    '[28,30,33,36,43,56;]'%2
    '[1;]'%4
    ''%5
    '[16;]'%6
    ''%7
    '[1,4,13,30,31,38,43,49,53,54;]'%8
    '[48;]'%9
    '[3;]'%10
    '[28;]'%11
    '[3,4,5,30,42;]'%12
    '[2,3,7,15,23,30,31,44,51,55;]'%13
    '[11,18;]'%14
    ''%15
    '[16,30;]'%16
    '[1,2,10,11,30,32,33,34,43,44;]'%17
    };
end
 
for sub = 1:length(subjects)
    subject = subjects(sub);
    subject = subject{:};
    

    load ([subject '_data_matrix_dim1']);
    d1 = x;
    load ([subject '_data_matrix_dim2']);
    d2 = x;
    clear x
    x = [d1(1) d1(2) d1(3)+d2(3)];
    
    events(:,1) = ones(x(3),1);  % which block one to six
    if b_order(sub) == 1
    events(:,2) = [(ones(1,d1(3))*1) (ones(1,d2(3))*2)];
    elseif b_order(sub) == 2
    events(:,2) = [(ones(1,d1(3))*2) (ones(1,d2(3))*1)];
    end    
    
    events_mat = events;

    if auto_exclude ==1
         fname= [subject '_total_data_ICA_ca'];
         load(fname);

        x=size(total_data_ICA);
        Nelectrodes = x(1);
        Nsamples = x(2);
        Nevents = x(3);
        total_data_sing = reshape(total_data_ICA,1,Nelectrodes*Nevents*Nsamples);
        SD_data = std(total_data_sing);
        
        total_data_trial = permute(total_data_ICA,[3 1 2]);
        total_data_trial = reshape(total_data_trial,Nevents,Nelectrodes*Nsamples);
        
        minexclude = [];
        for j = 1:Nevents
            if all(abs(total_data_trial(j,:)) < SD_data*SD_min)
                minexclude(length(minexclude)+1)=j;
            end
        end
        
        if ~isempty(minexclude)
            include = 1:Nevents;
            include(minexclude) = [];
            total_data_sing = reshape(total_data_ICA(:,:,include),1,Nelectrodes*(length(include))*Nsamples);
            SD_data = std(total_data_sing);
        end
        
        exclude = [];
        for j = 1:Nevents
            if any(abs(total_data_trial(j,:)) > SD_data*SD_thresh)
                exclude(length(exclude)+1)=j;
            end
        end

         excludes{sub,1} = sort(unique([exclude minexclude]),'ascend');
    end
    
    reject = [];
    if exclude_none ==1
        reject = '';
    else
        reject = excludes{sub,1};
        if ischar(reject)
            reject=str2num(reject);
        end
    end
    
 
events_mat(:,3)=ones(length(events_mat),1);
events_mat(reject,3)=zeros(length(reject),1);
 
ss= '_block_info';
ss = [char(subjects(sub)) ss]
save(ss, 'events_mat');

clear events

end
