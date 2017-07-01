subjects = {'H3', 'F1', 'F7', 'F6', 'OA2'};

for subject = 1:length(subjects)


d1 = dir([char(subjects(subject)) '_total_data_ica.mat']);
d2 = dir([char(subjects(subject)) ' _block_info.mat']);
load(d1.name); load(d2.name);

for i=1:size(total_data_ICA2,3)
    total_data_ICA2(:,:,i)=(detrend4(squeeze(total_data_ICA2(:,:,i))',1,125,125))';
end

total_data_ICA2 = blcorrect3(total_data_ICA2, 438);

block_size = max(events_mat(:,1));

for block = 1:block_size

    for trial = 1:6
    ii = find(events_mat(:,1)==block & events_mat(:,2)==trial & events_mat(:,3)==1);
avg = squeeze(mean(total_data_ICA2(:,:,ii),3));

name=[char(subjects(subject)) '_' num2str(trial) '_avgP2'];
%name=[char(subjects(subject)) '_' num2str(block) '_' num2str(trial) '_avg']; %%% use if there is more than one block
save(name, 'avg'); clear avg;

    end
end
end