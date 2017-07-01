%Peak definitions
% P2 - most positive peak between 200ms and 700ms
% N2 - most negative peak within 200ms prior to P2
% N1 - most negative peak between 150ms and N2 peak.

clear all

subjects = {'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13' 'H14', 'H15', 'H16', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'OA1', 'OA2', 'OA4', 'OA5', 'OA6', 'OA7', 'OA8', 'OA9', 'OA10', 'OA11', 'OA12','OA13','OA14', 'OA15', 'OA16','OA17'};
  
Nsub = length(subjects);
ncond = 6;
npeaks=5;
nchan = 64;
filt_spn_low = 0;
filt_spn_high = 20;
filt_lep_low = 2;
filt_lep_high = 20;

AMP = zeros(Nsub,ncond,npeaks,1); % subject, condition, time bins, electrode
LATEN = zeros(Nsub,npeaks,ncond);

elec_avg_early = [21 55 58]; % F1,Fz,F2
elec_avg_late = [14 15 20]; % Cz,CPz,FCz
%elec_avg_n1 = [5 39 41 47 25 27]; % C3,FC3,CP3,C5,FC5,CP5,
elec_avg_n1 = [18 24 25 47 53 57];% T3,FT3,TP3,C5,FC5,CP5,
elec_avg_n2 = [14 15 20]; % Fz,Cz,CPz
%elec_avg_p2 = [34 18 33]; % Cz,CPz,FCz
elec_avg_p2 = [14 15 20]; % Cz,CPz,FCz

for n =1:Nsub
    subject = subjects(n);
    subject = char(subject);   

    grandavg=zeros(nchan,2750);
    for c = 1:ncond
        fname= [subject '_' num2str(c) '_avg_ca'];
        eval(['mdata' num2str(c) '= load(fname);']);     
        eval(['mdata' num2str(c) '.avg = eegfilt(mdata' num2str(c) '.avg,500,filt_spn_low,filt_spn_high);']);  %filter the data
        eval(['grandavg=grandavg+mdata' num2str(c) '.avg;']);
    end
    grandavg= grandavg/ncond;


%Extract and save data

b1=mdata1.avg(:,250:500)';
b2=mdata2.avg(:,250:500)';
b3=mdata3.avg(:,250:500)';
b4=mdata4.avg(:,250:500)';
b5=mdata3.avg(:,250:500)';
b6=mdata4.avg(:,250:500)';
b1_av=mean(b1);
b2_av=mean(b2);
b3_av=mean(b3);
b4_av=mean(b4);
b5_av=mean(b5);
b6_av=mean(b6);
eval(['save ' subject '_1_b_av.dat'  ' b1_av /ASCII'])
eval(['save ' subject '_2_b_av.dat'  ' b2_av /ASCII'])
eval(['save ' subject '_3_b_av.dat'  ' b3_av /ASCII'])
eval(['save ' subject '_4_b_av.dat'  ' b4_av /ASCII'])
eval(['save ' subject '_5_b_av.dat'  ' b5_av /ASCII'])
eval(['save ' subject '_6_b_av.dat'  ' b6_av /ASCII'])


e1=mdata1.avg(:,750:1000)';
e2=mdata2.avg(:,750:1000)';
e3=mdata3.avg(:,750:1000)';
e4=mdata4.avg(:,750:1000)';
e5=mdata5.avg(:,750:1000)';
e6=mdata6.avg(:,750:1000)';
e1_av=mean(e1);
e2_av=mean(e2);
e3_av=mean(e3);
e4_av=mean(e4);
e5_av=mean(e5);
e6_av=mean(e6);
%c_e_av=[e4_av,e5_av,e6_av];
%c_e_av=mean(c_e_av,2);
%u_e_av=[e1_av,e2_av,e3_av];
%u_e_av=mean(u_e_av,2);
eval(['save ' subject '_1_e_av.dat'  ' e1_av /ASCII'])
eval(['save ' subject '_2_e_av.dat'  ' e2_av /ASCII'])
eval(['save ' subject '_3_e_av.dat'  ' e3_av /ASCII'])
eval(['save ' subject '_4_e_av.dat'  ' e4_av /ASCII'])
eval(['save ' subject '_5_e_av.dat'  ' e5_av /ASCII'])
eval(['save ' subject '_6_e_av.dat'  ' e6_av /ASCII'])
AMP(n,1,1,1) = mean(e1_av(:, elec_avg_early),2);
AMP(n,2,1,1) = mean(e2_av(:, elec_avg_early),2);
AMP(n,3,1,1) = mean(e3_av(:, elec_avg_early),2);
AMP(n,4,1,1) = mean(e4_av(:, elec_avg_early),2);
AMP(n,5,1,1) = mean(e5_av(:, elec_avg_early),2);
AMP(n,6,1,1) = mean(e6_av(:, elec_avg_early),2);


l1=mdata1.avg(:,1750:2000)';%late
l2=mdata2.avg(:,1750:2000)';
l3=mdata3.avg(:,1750:2000)';
l4=mdata4.avg(:,1750:2000)';
l5=mdata5.avg(:,1750:2000)';
l6=mdata6.avg(:,1750:2000)';
l1_av=mean(l1);
l2_av=mean(l2);
l3_av=mean(l3);
l4_av=mean(l4);
l5_av=mean(l5);
l6_av=mean(l6);
%c_l_av=[l4_av,l5_av,l6_av];
%c_l_av=mean(c_l_av,2);
%u_l_av=[l1_av,l2_av,l3_av];
%u_l_av=mean(u_l_av,2);

eval(['save ' subject '_1_l_av.dat'  ' l1_av /ASCII'])
eval(['save ' subject '_2_l_av.dat'  ' l2_av /ASCII'])
eval(['save ' subject '_3_l_av.dat'  ' l3_av /ASCII'])
eval(['save ' subject '_4_l_av.dat'  ' l4_av /ASCII'])
eval(['save ' subject '_5_l_av.dat'  ' l5_av /ASCII'])
eval(['save ' subject '_6_l_av.dat'  ' l6_av /ASCII'])
AMP(n,1,2,1) = mean(l1_av(:, elec_avg_late),2);
AMP(n,2,2,1) = mean(l2_av(:, elec_avg_late),2);
AMP(n,3,2,1) = mean(l3_av(:, elec_avg_late),2);
AMP(n,4,2,1) = mean(l4_av(:, elec_avg_late),2);
AMP(n,5,2,1) = mean(l5_av(:, elec_avg_late),2);
AMP(n,6,2,1) = mean(l6_av(:, elec_avg_late),2);

grandavg=zeros(nchan,2750);
for c = 1:ncond
    fname= [subject '_' num2str(c) '_avg_ca'];
    eval(['mdata' num2str(c) '= load(fname);']);     
    eval(['mdata' num2str(c) '.avg = eegfilt(mdata' num2str(c) '.avg,500,filt_lep_low,filt_lep_high);']);  %filter the data
    eval(['mdata' num2str(c) '.avg = blcorrect4(mdata' num2str(c) '.avg, 1750);']);
    eval(['grandavg=grandavg+mdata' num2str(c) '.avg;']);
end
grandavg= grandavg/ncond;


l1=mdata1.avg(:,1750:2000)';%late
l2=mdata2.avg(:,1750:2000)';
l3=mdata3.avg(:,1750:2000)';
l4=mdata4.avg(:,1750:2000)';
l5=mdata5.avg(:,1750:2000)';
l6=mdata6.avg(:,1750:2000)';
l1_av=mean(l1);
l2_av=mean(l2);
l3_av=mean(l3);
l4_av=mean(l4);
l5_av=mean(l5);
l6_av=mean(l6);
eval(['save ' subject '_1_b2_av.dat'  ' l1_av /ASCII'])
eval(['save ' subject '_2_b2_av.dat'  ' l2_av /ASCII'])
eval(['save ' subject '_3_b2_av.dat'  ' l3_av /ASCII'])
eval(['save ' subject '_4_b2_av.dat'  ' l4_av /ASCII'])
eval(['save ' subject '_5_b2_av.dat'  ' l5_av /ASCII'])
eval(['save ' subject '_6_b2_av.dat'  ' l6_av /ASCII'])


%identify latencies

[XMAX,IMAX,XMIN,IMIN] = extrema(mean(grandavg(elec_avg_p2,2100:2300),1));%%search for extrema on internet - peak detection prog, xmax ampli imax latency
LATp2 = IMAX(1)+2099; %% change to be the same as the 1st time point to zero everything
[XMAX,IMAX,XMIN,IMIN] = extrema(mean(grandavg(elec_avg_n2,LATp2-75:LATp2),1));%%search for extrema on internet - peak detection prog, xmax ampli imax latency
LATn2 = IMIN(1)+LATp2-76;

n1_range = max([2075 LATn2-100]):LATn2;
onesind = mean(-grandavg(elec_avg_n1,n1_range),1) > mean(-grandavg(elec_avg_n2,n1_range),1);
n1_grandavg = grandavg(:,n1_range).*repmat(onesind,size(grandavg,1),1);

[IMINa,XMINa] = findpeaks(mean(-n1_grandavg(elec_avg_n1,:),1));%%search for extrema on internet - peak detection prog, xmax ampli imax latency

if isempty(IMINa)
    [IMINa,XMINa] = findpeaks(mean(-grandavg(elec_avg_n1,n1_range),1));
end

top=0;
topi=find(XMINa==max(XMINa));%length(IMINa);

LATn1 = IMINa(topi)+max([2075 LATn2-100])-1;


[XMAX,IMAX,XMIN,IMIN] = extrema(mean(mdata1.avg(elec_avg_p2,LATp2-25:LATp2+25),1));%%each cond seperatly +/-64 ms 15 is CZ
LATp2_1 = IMAX(1)+LATp2-26;
[XMAX,IMAX,XMIN,IMIN] = extrema(mean(mdata2.avg(elec_avg_p2,LATp2-25:LATp2+25),1));
LATp2_2 = IMAX(1)+LATp2-26;
[XMAX,IMAX,XMIN,IMIN] = extrema(mean(mdata3.avg(elec_avg_p2,LATp2-25:LATp2+25),1));%%each cond seperatly +/-64 ms 15 is CZ
LATp2_3 = IMAX(1)+LATp2-26;
[XMAX,IMAX,XMIN,IMIN] = extrema(mean(mdata4.avg(elec_avg_p2,LATp2-25:LATp2+25),1));
LATp2_4 = IMAX(1)+LATp2-26;
[XMAX,IMAX,XMIN,IMIN] = extrema(mean(mdata5.avg(elec_avg_p2,LATp2-25:LATp2+25),1));%%each cond seperatly +/-64 ms 15 is CZ
LATp2_5 = IMAX(1)+LATp2-26;
[XMAX,IMAX,XMIN,IMIN] = extrema(mean(mdata6.avg(elec_avg_p2,LATp2-25:LATp2+25),1));
LATp2_6 = IMAX(1)+LATp2-26;

[XMAX,IMAX,XMIN,IMIN] = extrema(mean(mdata1.avg(elec_avg_n2,LATn2-25:LATn2+25),1));
LATn2_1 = IMIN(1)+LATn2-26;
[XMAX,IMAX,XMIN,IMIN] = extrema(mean(mdata2.avg(elec_avg_n2,LATn2-25:LATn2+25),1));
LATn2_2 = IMIN(1)+LATn2-26;
[XMAX,IMAX,XMIN,IMIN] = extrema(mean(mdata3.avg(elec_avg_n2,LATn2-25:LATn2+25),1));
LATn2_3 = IMIN(1)+LATn2-26;
[XMAX,IMAX,XMIN,IMIN] = extrema(mean(mdata4.avg(elec_avg_n2,LATn2-25:LATn2+25),1));
LATn2_4 = IMIN(1)+LATn2-26;
[XMAX,IMAX,XMIN,IMIN] = extrema(mean(mdata5.avg(elec_avg_n2,LATn2-25:LATn2+25),1));
LATn2_5 = IMIN(1)+LATn2-26;
[XMAX,IMAX,XMIN,IMIN] = extrema(mean(mdata6.avg(elec_avg_n2,LATn2-25:LATn2+25),1));
LATn2_6 = IMIN(1)+LATn2-26;

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
    eval(['LATn1_' num2str(c) ' = IMIN(topi)+max([2075 LATn1-50])-1;']);
end


peaks = {'n2', 'n1', 'p2'};
amp_pos = [4 3 5];
lat_pos = [2 1 3];
for p = 1:length(peaks)
    for c = 1:ncond;
        eval(['c' num2str(c) '=mdata' num2str(c) '.avg(:,(LAT' peaks{p} '_' num2str(c) ')-2:(LAT' peaks{p} '_' num2str(c) ')+2);']);
        eval(['c' num2str(c) '_av=mean(c' num2str(c) ',2);']);
        eval(['save ' subject '_' num2str(c) '_' peaks{p} '_av.dat'  ' c' num2str(c) '_av /ASCII;']);
        eval(['LATEN(n,' num2str(c) ',' num2str(lat_pos(p)) ') = ((LAT' peaks{p} '_' num2str(c) ')*2)-4000;']);
        eval(['AMP(n,' num2str(c) ',' num2str(amp_pos(p)) ',1) = mean(c' num2str(c) '_av(elec_avg_' peaks{p} ',:),1);']);
    end
end


end

AMP2 = reshape(AMP,Nsub,30); % reshapes: column order = cond(1->6) for each time (1->6), for each electrode (1->1).
save amplitudes.mat AMP2
save amplitudes.dat AMP2 /ASCII
P2LAT = squeeze(LATEN(:,:,3));
N2LAT = squeeze(LATEN(:,:,2));
N1LAT = squeeze(LATEN(:,:,1));
save P2latencies.mat P2LAT
save N2latencies.mat N2LAT
save N1latencies.mat N1LAT