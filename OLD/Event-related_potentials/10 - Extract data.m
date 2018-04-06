%Peak definitions
% P2 - most positive peak between 200ms and 700ms
% N2 - most negative peak within 200ms prior to P2
% N1 - most negative peak between 150ms and N2 peak.

clear all

subjects = {'P1_';'P2_';'P5_';'P6_';'P7_';'P8_';'P14_';'P16_';'P19_';'P20_';'P22_';'P23_';'P24_';'P25_';'P27_';'P30_';'P31_';'P32_';'P33_';'P35_';'S1_';'S2_';'S3_';'S5_';'S6_';'S8_';'S9_';'S10_';'S11_';'S18_';'S20_'};




Nsub = length(subjects);
ncond = 2;
npeaks = 5;
nchan = 62;
filt_spn_low = 0;
filt_spn_high = 20;
filt_lep_low = 2;
filt_lep_high = 20;

AMP = zeros(Nsub,ncond,npeaks,1); % subject, condition, time bins, electrode
LATEN = zeros(Nsub,npeaks,ncond);

ele = [5,6,17,18,19]; % electrodes to average when searching for peaks

elec_avg_early = [17 31 32]; % F1,Fz,F2
elec_avg_late = [18 19 61]; % Cz,CPz,Pz
%elec_avg_n1 = [5 39 41 47 25 27]; % C3,FC3,CP3,C5,FC5,CP5,
elec_avg_n1 = [13 53 55 47 25 27];% T3,FT3,TP3,C5,FC5,CP5,
elec_avg_n2 = [17 18 61]; % Fz,Cz,CPz
%elec_avg_p2 = [34 18 33]; % Cz,C1,C2
elec_avg_p2 = [17 18 61]; % Fz,Cz,CPz


for n =1:Nsub
    subject = subjects(n);
    subject = char(subject);   

    grandavg=zeros(nchan,2750);
    for c = 1:ncond
        fname= [subject 'avg_' num2str(c) '_ca'];
        eval(['mdata' num2str(c) '= load(fname);']);     
        eval(['mdata' num2str(c) '.avg = eegfilt(mdata' num2str(c) '.avg,500,filt_spn_low,filt_spn_high);']);  %filter the data
        eval(['grandavg=grandavg+mdata' num2str(c) '.avg;']);
    end
    grandavg= grandavg/ncond;


%Extract and save data

b1=mdata1.avg(:,250:500)';
b2=mdata2.avg(:,250:500)';
b1_av=mean(b1);
b2_av=mean(b2);
eval(['save ' subject '1_b_av.dat'  ' b1_av /ASCII'])
eval(['save ' subject '2_b_av.dat'  ' b2_av /ASCII'])

e1=mdata1.avg(:,750:1000)';
e2=mdata2.avg(:,750:1000)';
e1_av=mean(e1);
e2_av=mean(e2);
eval(['save ' subject '1_e_av.dat'  ' e1_av /ASCII'])
eval(['save ' subject '2_e_av.dat'  ' e2_av /ASCII'])
AMP(n,1,1,1) = mean(e1_av(:, elec_avg_early),2);
AMP(n,2,1,1) = mean(e2_av(:, elec_avg_early),2);

l1=mdata1.avg(:,1750:2000)';%late
l2=mdata2.avg(:,1750:2000)';
l1_av=mean(l1);
l2_av=mean(l2);
eval(['save ' subject '1_l_av.dat'  ' l1_av /ASCII'])
eval(['save ' subject '2_l_av.dat'  ' l2_av /ASCII'])
AMP(n,1,2,1) = mean(l1_av(:, elec_avg_late),2);
AMP(n,2,2,1) = mean(l2_av(:, elec_avg_late),2);

grandavg=zeros(nchan,2750);
for c = 1:ncond
    fname= [subject 'avg_' num2str(c) '_ca'];
    eval(['mdata' num2str(c) '= load(fname);']);     
    eval(['mdata' num2str(c) '.avg = eegfilt(mdata' num2str(c) '.avg,500,filt_lep_low,filt_lep_high);']);  %filter the data
    eval(['mdata' num2str(c) '.avg = blcorrect4(mdata' num2str(c) '.avg, 1750);']);
    eval(['grandavg=grandavg+mdata' num2str(c) '.avg;']);
end
grandavg= grandavg/ncond;


l1=mdata1.avg(:,1750:2000)';%late
l2=mdata2.avg(:,1750:2000)';
l1_av=mean(l1);
l2_av=mean(l2);
eval(['save ' subject '1_b2_av.dat'  ' l1_av /ASCII'])
eval(['save ' subject '2_b2_av.dat'  ' l2_av /ASCII'])


%identify latencies

[XMAX,IMAX,XMIN,IMIN] = extrema(mean(grandavg(elec_avg_p2,2100:2300),1));%%search for extrema on internet - peak detection prog, xmax ampli imax latency
LATp2 = IMAX(1)+2099; %% change to be the same as the 1st time point to zero everything
[XMAX,IMAX,XMIN,IMIN] = extrema(mean(grandavg(elec_avg_n2,LATp2-75:LATp2),1));%%search for extrema on internet - peak detection prog, xmax ampli imax latency
LATn2 = IMIN(1)+LATp2-76;

n1_range = max([2075 LATn2-100]):LATn2;
onesind = mean(-grandavg(elec_avg_n1,n1_range),1) > mean(-grandavg(elec_avg_n2,n1_range),1);
n1_grandavg = grandavg(:,n1_range).*repmat(onesind,size(grandavg,1),1);

[IMINa,XMINa] = findpeaks(mean(-n1_grandavg(elec_avg_n1,:),1));%%search for extrema on internet - peak detection prog, xmax ampli imax latency
%[XMAX,IMAX,XMINb,IMINb] = extrema(mean(grandavg(elec_avg_n2,LATn2-100:LATn2-10),1));%%search for extrema on internet - peak detection prog, xmax ampli imax latency
%if (XMINa(1)>XMINb(1))
%    LATn1 = IMINb(1)+LATn2-101;
%elseif (XMINa(1)<XMINb(1))

if isempty(IMINa)
    [IMINa,XMINa] = findpeaks(mean(-grandavg(elec_avg_n1,n1_range),1));
end

top=0;
topi=find(XMINa==max(XMINa));%length(IMINa);
%while top==0 && topi > 0
    LATn1 = XMINa(topi)+max([2075 LATn2-100])-1;
%    if mean(grandavg(elec_avg_n1,LATn1),1) < mean(grandavg(elec_avg_n1,LATn1+1),1) && mean(grandavg(elec_avg_n1,LATn1),1) < mean(grandavg(elec_avg_n1,LATn1-1),1)
%        top=1;
%        break;
%    else 
%        topi=topi-1;
%    end
%end
%end

% for those subjects that work best using IMIN(1) instead of max(IMIN):
%if strcmp(subject, 'H7') == 1 
%    LATn = IMIN(1)+2079;
%end



[XMAX,IMAX,XMIN,IMIN] = extrema(mean(mdata1.avg(elec_avg_p2,LATp2-25:LATp2+25),1));%%each cond seperatly +/-64 ms 15 is CZ
LATp2_1 = IMAX(1)+LATp2-26;
[XMAX,IMAX,XMIN,IMIN] = extrema(mean(mdata2.avg(elec_avg_p2,LATp2-25:LATp2+25),1));
LATp2_2 = IMAX(1)+LATp2-26;
[XMAX,IMAX,XMIN,IMIN] = extrema(mean(mdata1.avg(elec_avg_n2,LATn2-25:LATn2+25),1));
LATn2_1 = IMIN(1)+LATn2-26;
[XMAX,IMAX,XMIN,IMIN] = extrema(mean(mdata2.avg(elec_avg_n2,LATn2-25:LATn2+25),1));
LATn2_2 = IMIN(1)+LATn2-26;
%if (XMINa(1)>XMINb(1))
%[XMAX,IMAX,XMIN,IMIN] = extrema(mean(mdata1.avg(elec_avg_n2,LATn1-10:LATn1+10),1));
%LATn1_1 = IMIN(1)+LATn1-11;
%[XMAX,IMAX,XMIN,IMIN] = extrema(mean(mdata2.avg(elec_avg_n2,LATn1-10:LATn1+10),1));
%LATn1_2 = IMIN(1)+LATn1-11;
%elseif (XMINa(1)<XMINb(1))

for c = 1:ncond
    eval(['avg=mdata' num2str(c) '.avg;']);
    eval(['lat=LATn2_' num2str(c) ';']);
    n1_range = max([2075 LATn1-50]):min([lat LATn1+50]);
    onesind = mean(-avg(elec_avg_n1,n1_range),1) > mean(-avg(elec_avg_n2,n1_range),1);
    n1_mdata = avg(:,n1_range).*repmat(onesind,size(avg,1),1);
    [IMIN,XMIN] = findpeaks(mean(-n1_mdata(elec_avg_n1,:),1));%%search for extrema on internet - peak detection prog, xmax ampli imax latency
    if isempty(IMIN)
        [IMIN,XMIN] = findpeaks(mean(-avg(elec_avg_n1,n1_range),1));
        topi=find(XMIN==max(XMIN));%length(IMINa);
    end
    if isempty(IMIN)
        [XMAX,IMAX,XMIN,IMIN] = extrema(mean(avg(elec_avg_n1,n1_range),1));
        topi=find(XMIN==min(XMIN));%length(IMINa);
    else
        topi=find(XMIN==max(XMIN));%length(IMINa);
    end
    eval(['LATn1_' num2str(c) ' = XMIN(topi)+max([2075 LATn1-50])-1;']);
end

peaks = {'n2', 'n1', 'p2'};
amp_pos = [4 3 5];
lat_pos = [2 1 3];
for p = 1:length(peaks)
    for c = 1:ncond
        eval(['c' num2str(c) '=mdata' num2str(c) '.avg(:,(LAT' peaks{p} '_' num2str(c) ')-2:(LAT' peaks{p} '_' num2str(c) ')+2);']);
        eval(['c' num2str(c) '_av=mean(c' num2str(c) ',2);']);
        eval(['save ' subject '_' num2str(c) '_' peaks{p} '_av.dat'  ' c' num2str(c) '_av /ASCII;']);
        eval(['LATEN(n,' num2str(c) ',' num2str(lat_pos(p)) ') = ((LAT' peaks{p} '_' num2str(c) ')*2)-4000;']);
        eval(['AMP(n,' num2str(c) ',' num2str(amp_pos(p)) ',1) = mean(c' num2str(c) '_av(elec_avg_' peaks{p} ',:),1);']);
    end
end

end

AMP2 = reshape(AMP,Nsub,10); % reshapes: column order = cond(1->6) for each time (1->6), for each electrode (1->1).
save amplitudes.mat AMP2
save amplitudes.dat AMP2 /ASCII
P2LAT = squeeze(LATEN(:,:,3));
N2LAT = squeeze(LATEN(:,:,2));
N1LAT = squeeze(LATEN(:,:,1));
save P2latencies.mat P2LAT
save N2latencies.mat N2LAT
save N1latencies.mat N1LAT