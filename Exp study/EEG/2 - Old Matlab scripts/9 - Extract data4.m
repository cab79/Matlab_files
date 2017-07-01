
clear all

subjects = {'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13' 'H14', 'H15', 'H16', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'OA1', 'OA2', 'OA4', 'OA5', 'OA6', 'OA7', 'OA8', 'OA9', 'OA10', 'OA11', 'OA12','OA13','OA14', 'OA15', 'OA16','OA17'};
  
Nsub = length(subjects);

AMP = zeros(Nsub,6,6,1); % subject, condition, time bins, electrode
LATEN = zeros(Nsub,6,2);

elec_avg_early = [49 50 48 39 14 15 54 51 20];%[20 51 54];
elec_avg_mid = [49 50 48 39 14 15 54 51 20];%[15 20 49];
elec_avg_late = [49 50 48 39 14 15 54 51 20];%[14 15 49];
elec_avg_n2a = [49 50 48 39 14 15 54 51 20];%[14 15 49];
elec_avg_n2b = [49 50 48 39 14 15 54 51 20];%[13 19 47];
elec_avg_p2 = [49 50 48 39 14 15 54 51 20];%[15 20 51];

fnames={'_1_avg_ca.mat';
    '_2_avg_ca.mat';
    '_3_avg_ca.mat';
    '_4_avg_ca.mat';
    '_5_avg_ca.mat';
    '_6_avg_ca.mat';
  };

for n =1:Nsub
    subject = subjects(n);
    subject = char(subject);   

fname=char(fnames(1));
fname= [subject fname(1:9)];
mdata1= load(fname);       

fname=char(fnames(2));
fname= [subject fname(1:9)];
mdata2= load(fname);   

fname=char(fnames(3));
fname= [subject fname(1:9)];
mdata3= load(fname);   

fname=char(fnames(4));
fname= [subject fname(1:9)];
mdata4= load(fname);   

fname=char(fnames(5));
fname= [subject fname(1:9)];
mdata5= load(fname);   

fname=char(fnames(6));
fname= [subject fname(1:9)];
mdata6= load(fname);

    
grandavg= (mdata1.avg+mdata2.avg+mdata3.avg+mdata4.avg+mdata5.avg+mdata6.avg)/6;

Cz= grandavg(15, :);

%identify P2 latencies

[XMAX,IMAX,XMIN,IMIN] = extrema(Cz(:,2080:2360));%%search for extrema on internet - peak detection prog, xmax ampli imax latency
LAT = IMAX(1)+2079; %% change to be the same as the 1st time point to zero everything
[XMAX,IMAX,XMIN,IMIN] = extrema(Cz(:,2080:2200));%%search for extrema on internet - peak detection prog, xmax ampli imax latency
LATn = max(IMIN)+2079;

% for those subjects that work best using IMIN(1) instead of max(IMIN):
if strcmp(subject, 'H7') == 1 
    LATn = IMIN(1)+2079;
elseif strcmp(subject, 'H16') == 1 
    LATn = IMIN(1)+2079;
elseif strcmp(subject, 'F4') == 1 
    LATn = IMIN(1)+2079;
elseif strcmp(subject, 'F6') == 1 
    LATn = IMIN(1)+2079;
elseif strcmp(subject, 'F14') == 1 
    LATn = IMIN(1)+2079;
elseif strcmp(subject, 'F16') == 1 
    LATn = IMIN(1)+2079;
elseif strcmp(subject, 'OA7') == 1 
    LATn = IMIN(1)+2079;
end

[XMAX,IMAX,XMIN,IMIN] = extrema(mdata1.avg(15,LAT-60:LAT+60));%%each cond seperatly +/-64 ms 15 is CZ
LAT1 = IMAX(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata1.avg(15,LATn-60:LATn+60));
LATn1 = IMIN(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata2.avg(15,LAT-60:LAT+60));
LAT2 = IMAX(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata2.avg(15,LATn-60:LATn+60));
LATn2 = IMIN(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata3.avg(15,LAT-60:LAT+60));
LAT3 = IMAX(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata3.avg(15,LATn-60:LATn+60));
LATn3 = IMIN(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata4.avg(15,LAT-60:LAT+60));
LAT4 = IMAX(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata4.avg(15,LATn-60:LATn+60));
LATn4 = IMIN(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata5.avg(15,LAT-60:LAT+60));
LAT5 = IMAX(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata5.avg(15,LATn-60:LATn+60));
LATn5 = IMIN(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata6.avg(15,LAT-60:LAT+60));
LAT6 = IMAX(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata6.avg(15,LATn-60:LATn+60));
LATn6 = IMIN(1);

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




m1=mdata1.avg(:,1250:1500)';%mid time period
m2=mdata2.avg(:,1250:1500)';
m3=mdata3.avg(:,1250:1500)';
m4=mdata4.avg(:,1250:1500)';
m5=mdata5.avg(:,1250:1500)';
m6=mdata6.avg(:,1250:1500)';
m1_av=mean(m1);
m2_av=mean(m2);
m3_av=mean(m3);
m4_av=mean(m4);
m5_av=mean(m5);
m6_av=mean(m6);
%c_m_av=[m4_av,m5_av,m6_av];
%c_m_av=mean(c_m_av,2);
%u_m_av=[m1_av,m2_av,m3_av];
%u_m_av=mean(u_m_av,2);
eval(['save ' subject '_1_m_av.dat'  ' m1_av /ASCII'])
eval(['save ' subject '_2_m_av.dat'  ' m2_av /ASCII'])
eval(['save ' subject '_3_m_av.dat'  ' m3_av /ASCII'])
eval(['save ' subject '_4_m_av.dat'  ' m4_av /ASCII'])
eval(['save ' subject '_5_m_av.dat'  ' m5_av /ASCII'])
eval(['save ' subject '_6_m_av.dat'  ' m6_av /ASCII'])
AMP(n,1,2,1) = mean(m1_av(:, elec_avg_mid),2);
AMP(n,2,2,1) = mean(m2_av(:, elec_avg_mid),2);
AMP(n,3,2,1) = mean(m3_av(:, elec_avg_mid),2);
AMP(n,4,2,1) = mean(m4_av(:, elec_avg_mid),2);
AMP(n,5,2,1) = mean(m5_av(:, elec_avg_mid),2);
AMP(n,6,2,1) = mean(m6_av(:, elec_avg_mid),2);



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
AMP(n,1,3,1) = mean(l1_av(:, elec_avg_late),2);
AMP(n,2,3,1) = mean(l2_av(:, elec_avg_late),2);
AMP(n,3,3,1) = mean(l3_av(:, elec_avg_late),2);
AMP(n,4,3,1) = mean(l4_av(:, elec_avg_late),2);
AMP(n,5,3,1) = mean(l5_av(:, elec_avg_late),2);
AMP(n,6,3,1) = mean(l6_av(:, elec_avg_late),2);


mdata1.avg = blcorrect4(mdata1.avg, 1750);
mdata2.avg = blcorrect4(mdata2.avg, 1750);
mdata3.avg = blcorrect4(mdata3.avg, 1750);
mdata4.avg = blcorrect4(mdata4.avg, 1750);
mdata5.avg = blcorrect4(mdata5.avg, 1750);
mdata6.avg = blcorrect4(mdata6.avg, 1750);


n21=mdata1.avg(:,(LATn-61+LATn1)-2:(LATn-61+LATn1)+2)';%n2 for condtion 1 
n22=mdata2.avg(:,(LATn-61+LATn2)-2:(LATn-61+LATn2)+2)';%n2 for condition 2
n23=mdata3.avg(:,(LATn-61+LATn3)-2:(LATn-61+LATn3)+2)';
n24=mdata4.avg(:,(LATn-61+LATn4)-2:(LATn-61+LATn4)+2)';
n25=mdata5.avg(:,(LATn-61+LATn5)-2:(LATn-61+LATn5)+2)';
n26=mdata6.avg(:,(LATn-61+LATn6)-2:(LATn-61+LATn6)+2)';
n21_av=mean(n21);
n22_av=mean(n22);
n23_av=mean(n23);
n24_av=mean(n24);
n25_av=mean(n25);
n26_av=mean(n26);
%c_n2_av=[n24_av,n25_av,n26_av];
%c_n2_av=mean(c_n2_av,2);
%u_n2_av=[n21_av,n22_av,n23_av];
%u_n2_av=mean(u_n2_av,2);
eval(['save ' subject '_1_n2_av.dat'  ' n21_av /ASCII'])%all electrodes in 1 file - source analysis
eval(['save ' subject '_2_n2_av.dat'  ' n22_av /ASCII'])
eval(['save ' subject '_3_n2_av.dat'  ' n23_av /ASCII'])
eval(['save ' subject '_4_n2_av.dat'  ' n24_av /ASCII'])
eval(['save ' subject '_5_n2_av.dat'  ' n25_av /ASCII'])
eval(['save ' subject '_6_n2_av.dat'  ' n26_av /ASCII'])
LATEN(n,1,2) = ((LATn-61+LATn1)*2)-4000;
LATEN(n,2,2) = ((LATn-61+LATn2)*2)-4000;
LATEN(n,3,2) = ((LATn-61+LATn3)*2)-4000;
LATEN(n,4,2) = ((LATn-61+LATn4)*2)-4000;
LATEN(n,5,2) = ((LATn-61+LATn5)*2)-4000;
LATEN(n,6,2) = ((LATn-61+LATn6)*2)-4000;
AMP(n,1,4,1) = mean(n21_av(:, elec_avg_n2a),2);
AMP(n,2,4,1) = mean(n22_av(:, elec_avg_n2a),2);
AMP(n,3,4,1) = mean(n23_av(:, elec_avg_n2a),2);
AMP(n,4,4,1) = mean(n24_av(:, elec_avg_n2a),2);
AMP(n,5,4,1) = mean(n25_av(:, elec_avg_n2a),2);
AMP(n,6,4,1) = mean(n26_av(:, elec_avg_n2a),2);
AMP(n,1,5,1) = mean(n21_av(:, elec_avg_n2b),2);
AMP(n,2,5,1) = mean(n22_av(:, elec_avg_n2b),2);
AMP(n,3,5,1) = mean(n23_av(:, elec_avg_n2b),2);
AMP(n,4,5,1) = mean(n24_av(:, elec_avg_n2b),2);
AMP(n,5,5,1) = mean(n25_av(:, elec_avg_n2b),2);
AMP(n,6,5,1) = mean(n26_av(:, elec_avg_n2b),2);

p21=mdata1.avg(:,(LAT-61+LAT1)-2:(LAT-61+LAT1)+2)';%P2 for condtion 1 
p22=mdata2.avg(:,(LAT-61+LAT2)-2:(LAT-61+LAT2)+2)';%P2 for condition 2
p23=mdata3.avg(:,(LAT-61+LAT3)-2:(LAT-61+LAT3)+2)';
p24=mdata4.avg(:,(LAT-61+LAT4)-2:(LAT-61+LAT4)+2)';
p25=mdata5.avg(:,(LAT-61+LAT5)-2:(LAT-61+LAT5)+2)';
p26=mdata6.avg(:,(LAT-61+LAT6)-2:(LAT-61+LAT6)+2)';
p21_av=mean(p21);
p22_av=mean(p22);
p23_av=mean(p23);
p24_av=mean(p24);
p25_av=mean(p25);
p26_av=mean(p26);
AMP(n,1,6,1) = mean(p21_av(:, elec_avg_p2),2);
AMP(n,2,6,1) = mean(p22_av(:, elec_avg_p2),2);
AMP(n,3,6,1) = mean(p23_av(:, elec_avg_p2),2);
AMP(n,4,6,1) = mean(p24_av(:, elec_avg_p2),2);
AMP(n,5,6,1) = mean(p25_av(:, elec_avg_p2),2);
AMP(n,6,6,1) = mean(p26_av(:, elec_avg_p2),2);
eval(['save ' subject '_1_p2_av.dat'  ' p21_av /ASCII'])%all electrodes in 1 file - source analysis
eval(['save ' subject '_2_p2_av.dat'  ' p22_av /ASCII'])
eval(['save ' subject '_3_p2_av.dat'  ' p23_av /ASCII'])
eval(['save ' subject '_4_p2_av.dat'  ' p24_av /ASCII'])
eval(['save ' subject '_5_p2_av.dat'  ' p25_av /ASCII'])
eval(['save ' subject '_6_p2_av.dat'  ' p26_av /ASCII'])
LATEN(n,1,1) = ((LAT-61+LAT1)*2)-4000;
LATEN(n,2,1) = ((LAT-61+LAT2)*2)-4000;
LATEN(n,3,1) = ((LAT-61+LAT3)*2)-4000;
LATEN(n,4,1) = ((LAT-61+LAT4)*2)-4000;
LATEN(n,5,1) = ((LAT-61+LAT5)*2)-4000;
LATEN(n,6,1) = ((LAT-61+LAT6)*2)-4000;

end

AMP2 = reshape(AMP,Nsub,36); % reshapes: column order = cond(1->6) for each time (1->6), for each electrode (1->1).
save amplitudes.mat AMP2
save amplitudes.dat AMP2 /ASCII
P2LAT = squeeze(LATEN(:,:,1));
N2LAT = squeeze(LATEN(:,:,2));
save P2latencies.mat P2LAT
save N2latencies.mat N2LAT