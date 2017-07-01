clear all

subjects = {'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13' 'H14', 'H15', 'H16', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'OA1', 'OA2', 'OA4', 'OA5', 'OA6', 'OA7', 'OA8', 'OA9', 'OA10', 'OA11', 'OA12', 'OA13', 'OA14', 'OA15', 'OA16', 'OA17'};


gamma1 = zeros(64,length(subjects),2);
gamma2 = zeros(64,length(subjects),2);
gamma3 = zeros(64,length(subjects),2);
theta = zeros(64,length(subjects),2);
alpha1 = zeros(64,length(subjects),2);
alpha2 = zeros(64,length(subjects),2);
alpha3 = zeros(64,length(subjects),2);
beta = zeros(64,length(subjects),2);
delta = zeros(64,length(subjects),2);

IAFs=[];

for subject = 1:length(subjects)

    


d1 = [char(subjects(subject)) '_total_data_ICA_ca.mat'];
d2 = [char(subjects(subject)) '_block_info.mat'];
load(d1); load(d2); 

ii1 = find(events_mat(:,2) == 1 & events_mat(:,3)==1);
total_data_ICA1 = total_data_ICA(:,:,ii1);

ii2 = find(events_mat(:,2) == 2 & events_mat(:,3)==1);
total_data_ICA2 = total_data_ICA(:,:,ii2);

names = {'total_data_ICA1'
    'total_data_ICA2'
    };

for cond = 1:2
    iaf=[];
    eval(['data = ' char(names(cond)) ';' ]);

frames = 999 %frames per epoch {0 -> data length}
srate  = 500 %sampling rate per channel (Hz)

cond
[spectra,freqs,speccomp,contrib,specstd] = spectopo(data, frames, srate,'freqrange',[2 20],'freqfac',10);

ca = mean(spectra,1);
iaf = freqs(find(ca==max(ca(:,find(freqs>8 & freqs<12)))));
IAFs(subject,cond) = iaf;
IAFs(subject,cond+2) = ca(find(ca==max(ca(:,find(freqs>8 & freqs<12)))));

eval(['cond' num2str(cond) '.spectra = spectra;']);
eval(['cond' num2str(cond) '.freqs = freqs;']);

end

if IAFs(subject,3)>IAFs(subject,4); 
    temp1 = cond2;
    temp2 = cond1;
    cond1 = temp1;
    cond2 = temp2;
    iaftemp1 = IAFs(subject,2);
    iaftemp2 = IAFs(subject,1);
    IAFs(subject,1) = iaftemp1;
    IAFs(subject,2) = iaftemp2;
    iaftemp3 = IAFs(subject,4);
    iaftemp4 = IAFs(subject,3);
    IAFs(subject,3) = iaftemp3;
    IAFs(subject,4) = iaftemp4;
    events_mat(ii1,2) = 2*ones(length(ii1),1);
    events_mat(ii2,2) = ones(length(ii2),1);
end

for cond = 1:2
    
    eval(['spectra = cond' num2str(cond) '.spectra;']);
    eval(['freqs = cond' num2str(cond) '.freqs;']);
    
    g1 = find(freqs>30 & freqs<60);
    gamma1(:,subject,cond) = mean(spectra(:,g1),2);
    g2 = find(freqs>60 & freqs<90);
    gamma2(:,subject,cond) = mean(spectra(:,g2),2);
    g3 = find(freqs>90 & freqs<120);
    gamma3(:,subject,cond) = mean(spectra(:,g3),2);
    b = find(freqs>iaf+2 & freqs<20);
    beta(:,subject,cond) = mean(spectra(:,b),2);
    a1 = find(freqs>iaf-4 & freqs<iaf-2);
    alpha1(:,subject,cond) = mean(spectra(:,a1),2);
    a2 = find(freqs>iaf-2 & freqs<iaf);
    alpha2(:,subject,cond) = mean(spectra(:,a2),2);
    a3 = find(freqs>iaf & freqs<iaf+2);
    alpha3(:,subject,cond) = mean(spectra(:,a3),2);
    t = find(freqs>iaf-6 & freqs<iaf-4);
    theta(:,subject,cond) = mean(spectra(:,t),2);
    d = find(freqs>2 & freqs<4);
    delta(:,subject,cond) = mean(spectra(:,d),2);

end

ss= '_block_info2';
ss = [char(subjects(subject)) ss]
save(ss, 'events_mat');

IAFs

end

save gamma1.mat gamma1
save gamma2.mat gamma2
save gamma3.mat gamma3
save beta.mat beta
save alpha1.mat alpha1
save alpha2.mat alpha2
save alpha3.mat alpha3
save theta.mat theta
save delta.mat delta
save IAFs.mat IAFs