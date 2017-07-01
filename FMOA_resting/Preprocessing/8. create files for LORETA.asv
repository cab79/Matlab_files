clear all

subjects = {'M1';'M2';'M4';'M6';'M7';'M9';'M10';'M12';'M14';'M15';'M16';'M17';'M19';'M20';'M21';'M22';'M24';'M25';'M26';'M28';'M29';'M30';'M32';'M33';'M35';'M36';'M37';'M38';'M40'};


for subject = 1:length(subjects)

d1 = [char(subjects(subject)) '_total_data_ICAA.mat'];
d2 = [char(subjects(subject)) '_block_infoA.mat'];
load(d1); load(d2); 

ii1 = find(events_mat(:,2) == 1 & events_mat(:,3)==1);
total_data_ICAA1 = total_data_ICAA(:,:,ii1);

ii2 = find(events_mat(:,2) == 2 & events_mat(:,3)==1);
total_data_ICAA2 = total_data_ICAA(:,:,ii2);

d3 = [char(subjects(subject)) '_total_data_ICAB.mat'];
d4 = [char(subjects(subject)) '_block_infoB.mat'];
load(d3); load(d4); 

ii3 = find(events_mat(:,2) == 1 & events_mat(:,3)==1);
total_data_ICAB1 = total_data_ICAB(:,:,ii3);

ii4 = find(events_mat(:,2) == 2 & events_mat(:,3)==1);
total_data_ICAB2 = total_data_ICAB(:,:,ii4);

total_data_ICAA1 = reshape(total_data_ICAA1,size(total_data_ICAA1,1),size(total_data_ICAA1,2)*size(total_data_ICAA1,3));
total_data_ICAA2 = reshape(total_data_ICAA2,size(total_data_ICAA2,1),size(total_data_ICAA2,2)*size(total_data_ICAA2,3));
total_data_ICAB1 = reshape(total_data_ICAB1,size(total_data_ICAB1,1),size(total_data_ICAB1,2)*size(total_data_ICAB1,3));
total_data_ICAB2 = reshape(total_data_ICAB2,size(total_data_ICAB2,1),size(total_data_ICAB2,2)*size(total_data_ICAB2,3));

eval(['save ' num2str(subject) '_total_data_ICAA1.dat'  ' total_data_ICAA1 /ASCII']);
eval(['save ' num2str(subject) '_total_data_ICAA2.dat'  ' total_data_ICAA2 /ASCII']);
eval(['save ' num2str(subject) '_total_data_ICAB1.dat'  ' total_data_ICAB1 /ASCII']);
eval(['save ' num2str(subject) '_total_data_ICAB2.dat'  ' total_data_ICAB2 /ASCII']);

end