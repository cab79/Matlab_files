clear all
subjects = {'H3', 'F1','F4', 'F7', 'OA2'};

for subject = 1:length(subjects)


d1 = dir([char(subjects(subject)) '_total_data_ica.mat']);
d2 = dir([char(subjects(subject)) ' _block_info.mat']);
load(d1.name); load(d2.name);

for i=1:size(total_data_ICA2,3)
    total_data_ICA2(:,:,i)=(detrend2(squeeze(total_data_ICA2(:,:,i))',250,250))';
end

total_data_ICA2 = blcorrect4(total_data_ICA2, 250);

block_size = max(events_mat(:,1));

for block = 1:block_size

    for trial = 1:6
    ii = find(events_mat(:,1)==block & events_mat(:,2)==trial & events_mat(:,3)==1);
avg = squeeze(mean(total_data_ICA2(:,:,ii),3));

name=[char(subjects(subject)) '_' num2str(trial) '_avg'];
%name=[char(subjects(subject)) '_' num2str(block) '_' num2str(trial) '_avg']; %%% use if there is more than one block
save(name, 'avg'); clear avg;

    end
end
end




