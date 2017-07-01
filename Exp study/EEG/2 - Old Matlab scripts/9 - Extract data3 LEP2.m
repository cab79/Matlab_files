%extract data from C3, C5, CP3, CP5, average of (C3, C5, CP3, CP5)


clear all

subjects = {'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13' 'H14', 'H15', 'H16', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'OA1', 'OA2', 'OA4', 'OA5', 'OA6', 'OA7', 'OA8', 'OA9', 'OA10', 'OA11', 'OA12'};
  
Nsub = length(subjects);

AMP = zeros(Nsub,6,5,5); % subject, condition, time bins, electrode
LATEN = zeros(Nsub,6,2);

fnames={'_1_avgP2_ca.mat';
    '_2_avgP2_ca.mat';
    '_3_avgP2_ca.mat';
    '_4_avgP2_ca.mat';
    '_5_avgP2_ca.mat';
    '_6_avgP2_ca.mat';
  };

for n = 1:Nsub
    subject = subjects(n);
    subject = char(subject);   

fname=char(fnames(1));
fname= [subject fname(1:11)];
mdata1= load(fname);       
    
fname=char(fnames(2));
fname= [subject fname(1:11)];
mdata2= load(fname);   

fname=char(fnames(3));
fname= [subject fname(1:11)];
mdata3= load(fname);   

fname=char(fnames(4));
fname= [subject fname(1:11)];
mdata4= load(fname);   

fname=char(fnames(5));
fname= [subject fname(1:11)];
mdata5= load(fname);   

fname=char(fnames(6));
fname= [subject fname(1:11)];
mdata6= load(fname);

    
grandavg= (mdata1.avg+mdata2.avg+mdata3.avg+mdata4.avg+mdata5.avg+mdata6.avg)/6;

C3= grandavg(19, :);

%identify P2 latencies

[XMAX,IMAX,XMIN,IMIN] = extrema(C3(:,530:590));%%search for extrema on internet - peak detection prog, xmax ampli imax latency
LAT = IMAX(1)+529; %% change to be the same as the 1st time point to zero everything
[XMAX,IMAX,XMIN,IMIN] = extrema(C3(:,520:575));%%search for extrema on internet - peak detection prog, xmax ampli imax latency
LATn = IMIN(1)+519;

% for those subjects that work best using IMIN(1) instead of max(IMIN):
%if strcmp(subject, 'H7') == 1 
%    LATn = IMIN(1)+519;
%elseif strcmp(subject, 'H16') == 1 
%    LATn = IMIN(1)+519;
%elseif strcmp(subject, 'F4') == 1 
%    LATn = IMIN(1)+519;
%elseif strcmp(subject, 'F6') == 1 
%    LATn = IMIN(1)+519;
%elseif strcmp(subject, 'F14') == 1 
%    LATn = IMIN(1)+519;
%elseif strcmp(subject, 'F16') == 1 
%    LATn = IMIN(1)+519;
%elseif strcmp(subject, 'OA4') == 1 
%    LATn = IMIN(1)+519;
%elseif strcmp(subject, 'OA5') == 1 
%    LATn = IMIN(1)+519;
%elseif strcmp(subject, 'OA7') == 1 
%    [XMAX,IMAX,XMIN,IMIN] = extrema(C3(:,540:LAT)); %shifted the window forward a bit for this subject
%    LATn = IMIN(1)+539;
%elseif strcmp(subject, 'OA9') == 1 
%    LATn = IMIN(1)+519;
%elseif strcmp(subject, 'OA10') == 1 
%    LATn = IMIN(1)+519;
%end

[XMAX,IMAX,XMIN,IMIN] = extrema(mdata1.avg(19,LAT-8:LAT+8));%%each cond seperatly +/-64 ms 15 is CZ
LAT1 = IMAX(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata1.avg(19,LATn-8:LATn+8));
LATn1 = IMIN(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata2.avg(19,LAT-8:LAT+8));
LAT2 = IMAX(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata2.avg(19,LATn-8:LATn+8));
LATn2 = IMIN(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata3.avg(19,LAT-8:LAT+8));
LAT3 = IMAX(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata3.avg(19,LATn-8:LATn+8));
LATn3 = IMIN(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata4.avg(19,LAT-8:LAT+8));
LAT4 = IMAX(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata4.avg(19,LATn-8:LATn+8));
LATn4 = IMIN(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata5.avg(19,LAT-8:LAT+8));
LAT5 = IMAX(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata5.avg(19,LATn-8:LATn+8));
LATn5 = IMIN(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata6.avg(19,LAT-8:LAT+8));
LAT6 = IMAX(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata6.avg(19,LATn-8:LATn+8));
LATn6 = IMIN(1);

%Extract and save data

e1=mdata1.avg(:,188:250)';
e2=mdata2.avg(:,188:250)';
e3=mdata3.avg(:,188:250)';
e4=mdata4.avg(:,188:250)';
e5=mdata5.avg(:,188:250)';
e6=mdata6.avg(:,188:250)';
e1_av=mean(e1);
e2_av=mean(e2);
e3_av=mean(e3);
e4_av=mean(e4);
e5_av=mean(e5);
e6_av=mean(e6);
%c_e_av=[m4_av,m5_av,m6_av];
%c_e_av=mean(c_e_av,2);
%u_e_av=[m1_av,m2_av,m3_av];
%u_e_av=mean(u_e_av,2);
eval(['save ' subject '_1_e_av.dat'  ' e1_av /ASCII'])
eval(['save ' subject '_2_e_av.dat'  ' e2_av /ASCII'])
eval(['save ' subject '_3_e_av.dat'  ' e3_av /ASCII'])
eval(['save ' subject '_4_e_av.dat'  ' e4_av /ASCII'])
eval(['save ' subject '_5_e_av.dat'  ' e5_av /ASCII'])
eval(['save ' subject '_6_e_av.dat'  ' e6_av /ASCII'])


%save c_e_av.dat c_e_av /ASCII
%save u_e_av.dat u_e_av /ASCII
AMP(n,1,1,1) = e1_av(:, 19);
AMP(n,2,1,1) = e2_av(:, 19);
AMP(n,3,1,1) = e3_av(:, 19);
AMP(n,4,1,1) = e4_av(:, 19);
AMP(n,5,1,1) = e5_av(:, 19);
AMP(n,6,1,1) = e6_av(:, 19);
AMP(n,1,1,2) = e1_av(:, 53);
AMP(n,2,1,2) = e2_av(:, 53);
AMP(n,3,1,2) = e3_av(:, 53);
AMP(n,4,1,2) = e4_av(:, 53);
AMP(n,5,1,2) = e5_av(:, 53);
AMP(n,6,1,2) = e6_av(:, 53);
AMP(n,1,1,3) = e1_av(:, 13);
AMP(n,2,1,3) = e2_av(:, 13);
AMP(n,3,1,3) = e3_av(:, 13);
AMP(n,4,1,3) = e4_av(:, 13);
AMP(n,5,1,3) = e5_av(:, 13);
AMP(n,6,1,3) = e6_av(:, 13);
AMP(n,1,1,4) = e1_av(:, 47);
AMP(n,2,1,4) = e2_av(:, 47);
AMP(n,3,1,4) = e3_av(:, 47);
AMP(n,4,1,4) = e4_av(:, 47);
AMP(n,5,1,4) = e5_av(:, 47);
AMP(n,6,1,4) = e6_av(:, 47);
AMP(n,1,1,5) = mean([e1_av(:, 19),e1_av(:, 53),e1_av(:, 13),e1_av(:, 47)],2);
AMP(n,2,1,5) = mean([e2_av(:, 19),e2_av(:, 53),e2_av(:, 13),e2_av(:, 47)],2);
AMP(n,3,1,5) = mean([e3_av(:, 19),e3_av(:, 53),e3_av(:, 13),e3_av(:, 47)],2);
AMP(n,4,1,5) = mean([e4_av(:, 19),e4_av(:, 53),e4_av(:, 13),e4_av(:, 47)],2);
AMP(n,5,1,5) = mean([e5_av(:, 19),e5_av(:, 53),e5_av(:, 13),e5_av(:, 47)],2);
AMP(n,6,1,5) = mean([e6_av(:, 19),e6_av(:, 53),e6_av(:, 13),e6_av(:, 47)],2);

m1=mdata1.avg(:,313:375)';%mid time period
m2=mdata2.avg(:,313:375)';
m3=mdata3.avg(:,313:375)';
m4=mdata4.avg(:,313:375)';
m5=mdata5.avg(:,313:375)';
m6=mdata6.avg(:,313:375)';
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


AMP(n,1,2,1) = m1_av(:, 19);
AMP(n,2,2,1) = m2_av(:, 19);
AMP(n,3,2,1) = m3_av(:, 19);
AMP(n,4,2,1) = m4_av(:, 19);
AMP(n,5,2,1) = m5_av(:, 19);
AMP(n,6,2,1) = m6_av(:, 19);
AMP(n,1,2,2) = m1_av(:, 53);
AMP(n,2,2,2) = m2_av(:, 53);
AMP(n,3,2,2) = m3_av(:, 53);
AMP(n,4,2,2) = m4_av(:, 53);
AMP(n,5,2,2) = m5_av(:, 53);
AMP(n,6,2,2) = m6_av(:, 53);
AMP(n,1,2,3) = m1_av(:, 13);
AMP(n,2,2,3) = m2_av(:, 13);
AMP(n,3,2,3) = m3_av(:, 13);
AMP(n,4,2,3) = m4_av(:, 13);
AMP(n,5,2,3) = m5_av(:, 13);
AMP(n,6,2,3) = m6_av(:, 13);
AMP(n,1,2,4) = m1_av(:, 47);
AMP(n,2,2,4) = m2_av(:, 47);
AMP(n,3,2,4) = m3_av(:, 47);
AMP(n,4,2,4) = m4_av(:, 47);
AMP(n,5,2,4) = m5_av(:, 47);
AMP(n,6,2,4) = m6_av(:, 47);
AMP(n,1,2,5) = mean([m1_av(:, 19),m1_av(:, 53),m1_av(:, 13),m1_av(:, 47)],2);
AMP(n,2,2,5) = mean([m2_av(:, 19),m2_av(:, 53),m2_av(:, 13),m2_av(:, 47)],2);
AMP(n,3,2,5) = mean([m3_av(:, 19),m3_av(:, 53),m3_av(:, 13),m3_av(:, 47)],2);
AMP(n,4,2,5) = mean([m4_av(:, 19),m4_av(:, 53),m4_av(:, 13),m4_av(:, 47)],2);
AMP(n,5,2,5) = mean([m5_av(:, 19),m5_av(:, 53),m5_av(:, 13),m5_av(:, 47)],2);
AMP(n,6,2,5) = mean([m6_av(:, 19),m6_av(:, 53),m6_av(:, 13),m6_av(:, 47)],2);

l1=mdata1.avg(:,438:500)';%late
l2=mdata2.avg(:,438:500)';
l3=mdata3.avg(:,438:500)';
l4=mdata4.avg(:,438:500)';
l5=mdata5.avg(:,438:500)';
l6=mdata6.avg(:,438:500)';
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


%save c_l_av.dat c_l_av /ASCII
%save u_l_av.dat u_l_av /ASCII
AMP(n,1,3,1) = l1_av(:, 19);
AMP(n,2,3,1) = l2_av(:, 19);
AMP(n,3,3,1) = l3_av(:, 19);
AMP(n,4,3,1) = l4_av(:, 19);
AMP(n,5,3,1) = l5_av(:, 19);
AMP(n,6,3,1) = l6_av(:, 19);
AMP(n,1,3,2) = l1_av(:, 53);
AMP(n,2,3,2) = l2_av(:, 53);
AMP(n,3,3,2) = l3_av(:, 53);
AMP(n,4,3,2) = l4_av(:, 53);
AMP(n,5,3,2) = l5_av(:, 53);
AMP(n,6,3,2) = l6_av(:, 53);
AMP(n,1,3,3) = l1_av(:, 13);
AMP(n,2,3,3) = l2_av(:, 13);
AMP(n,3,3,3) = l3_av(:, 13);
AMP(n,4,3,3) = l4_av(:, 13);
AMP(n,5,3,3) = l5_av(:, 13);
AMP(n,6,3,3) = l6_av(:, 13);
AMP(n,1,3,4) = l1_av(:, 47);
AMP(n,2,3,4) = l2_av(:, 47);
AMP(n,3,3,4) = l3_av(:, 47);
AMP(n,4,3,4) = l4_av(:, 47);
AMP(n,5,3,4) = l5_av(:, 47);
AMP(n,6,3,4) = l6_av(:, 47);
AMP(n,1,3,5) = mean([l1_av(:, 19),l1_av(:, 53),l1_av(:, 13),l1_av(:, 47)],2);
AMP(n,2,3,5) = mean([l2_av(:, 19),l2_av(:, 53),l2_av(:, 13),l2_av(:, 47)],2);
AMP(n,3,3,5) = mean([l3_av(:, 19),l3_av(:, 53),l3_av(:, 13),l3_av(:, 47)],2);
AMP(n,4,3,5) = mean([l4_av(:, 19),l4_av(:, 53),l4_av(:, 13),l4_av(:, 47)],2);
AMP(n,5,3,5) = mean([l5_av(:, 19),l5_av(:, 53),l5_av(:, 13),l5_av(:, 47)],2);
AMP(n,6,3,5) = mean([l6_av(:, 19),l6_av(:, 53),l6_av(:, 13),l6_av(:, 47)],2);

p21=mdata1.avg(:,(LAT-9+LAT1)-1:(LAT-9+LAT1)+1)';%P2 for condtion 1 
p22=mdata2.avg(:,(LAT-9+LAT2)-1:(LAT-9+LAT2)+1)';%P2 for condition 2
p23=mdata3.avg(:,(LAT-9+LAT3)-1:(LAT-9+LAT3)+1)';
p24=mdata4.avg(:,(LAT-9+LAT4)-1:(LAT-9+LAT4)+1)';
p25=mdata5.avg(:,(LAT-9+LAT5)-1:(LAT-9+LAT5)+1)';
p26=mdata6.avg(:,(LAT-9+LAT6)-1:(LAT-9+LAT6)+1)';
p21_av=mean(p21);
p22_av=mean(p22);
p23_av=mean(p23);
p24_av=mean(p24);
p25_av=mean(p25);
p26_av=mean(p26);
%c_p2_av=[p24_av,p25_av,p26_av];
%c_p2_av=mean(c_p2_av,2);
%u_p2_av=[p21_av,p22_av,p23_av];
%u_p2_av=mean(u_p2_av,2);

eval(['save ' subject '_1_p2_av.dat'  ' p21_av /ASCII'])%all electrodes in 1 file - source analysis
eval(['save ' subject '_2_p2_av.dat'  ' p22_av /ASCII'])
eval(['save ' subject '_3_p2_av.dat'  ' p23_av /ASCII'])
eval(['save ' subject '_4_p2_av.dat'  ' p24_av /ASCII'])
eval(['save ' subject '_5_p2_av.dat'  ' p25_av /ASCII'])
eval(['save ' subject '_6_p2_av.dat'  ' p26_av /ASCII'])

LATEN(n,1,1) = ((LAT-9+LAT1)*8)-4000;
LATEN(n,2,1) = ((LAT-9+LAT2)*8)-4000;
LATEN(n,3,1) = ((LAT-9+LAT3)*8)-4000;
LATEN(n,4,1) = ((LAT-9+LAT4)*8)-4000;
LATEN(n,5,1) = ((LAT-9+LAT5)*8)-4000;
LATEN(n,6,1) = ((LAT-9+LAT6)*8)-4000;

AMP(n,1,4,1) = p21_av(:, 19);
AMP(n,2,4,1) = p22_av(:, 19);
AMP(n,3,4,1) = p23_av(:, 19);
AMP(n,4,4,1) = p24_av(:, 19);
AMP(n,5,4,1) = p25_av(:, 19);
AMP(n,6,4,1) = p26_av(:, 19);
AMP(n,1,4,2) = p21_av(:, 53);
AMP(n,2,4,2) = p22_av(:, 53);
AMP(n,3,4,2) = p23_av(:, 53);
AMP(n,4,4,2) = p24_av(:, 53);
AMP(n,5,4,2) = p25_av(:, 53);
AMP(n,6,4,2) = p26_av(:, 53);
AMP(n,1,4,3) = p21_av(:, 13);
AMP(n,2,4,3) = p22_av(:, 13);
AMP(n,3,4,3) = p23_av(:, 13);
AMP(n,4,4,3) = p24_av(:, 13);
AMP(n,5,4,3) = p25_av(:, 13);
AMP(n,6,4,3) = p26_av(:, 13);
AMP(n,1,4,4) = p21_av(:, 47);
AMP(n,2,4,4) = p22_av(:, 47);
AMP(n,3,4,4) = p23_av(:, 47);
AMP(n,4,4,4) = p24_av(:, 47);
AMP(n,5,4,4) = p25_av(:, 47);
AMP(n,6,4,4) = p26_av(:, 47);
AMP(n,1,4,5) = mean([p21_av(:, 19),p21_av(:, 53),p21_av(:, 13),p21_av(:, 47)],2);
AMP(n,2,4,5) = mean([p22_av(:, 19),p22_av(:, 53),p22_av(:, 13),p22_av(:, 47)],2);
AMP(n,3,4,5) = mean([p23_av(:, 19),p23_av(:, 53),p23_av(:, 13),p23_av(:, 47)],2);
AMP(n,4,4,5) = mean([p24_av(:, 19),p24_av(:, 53),p24_av(:, 13),p24_av(:, 47)],2);
AMP(n,5,4,5) = mean([p25_av(:, 19),p25_av(:, 53),p25_av(:, 13),p25_av(:, 47)],2);
AMP(n,6,4,5) = mean([p26_av(:, 19),p26_av(:, 53),p26_av(:, 13),p26_av(:, 47)],2);

n21=mdata1.avg(:,(LATn-9+LATn1)-1:(LATn-9+LATn1)+1)';%n2 for condtion 1 
n22=mdata2.avg(:,(LATn-9+LATn2)-1:(LATn-9+LATn2)+1)';%n2 for condition 2
n23=mdata3.avg(:,(LATn-9+LATn3)-1:(LATn-9+LATn3)+1)';
n24=mdata4.avg(:,(LATn-9+LATn4)-1:(LATn-9+LATn4)+1)';
n25=mdata5.avg(:,(LATn-9+LATn5)-1:(LATn-9+LATn5)+1)';
n26=mdata6.avg(:,(LATn-9+LATn6)-1:(LATn-9+LATn6)+1)';
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

LATEN(n,1,2) = ((LATn-9+LATn1)*8)-4000;
LATEN(n,2,2) = ((LATn-9+LATn2)*8)-4000;
LATEN(n,3,2) = ((LATn-9+LATn3)*8)-4000;
LATEN(n,4,2) = ((LATn-9+LATn4)*8)-4000;
LATEN(n,5,2) = ((LATn-9+LATn5)*8)-4000;
LATEN(n,6,2) = ((LATn-9+LATn6)*8)-4000;

AMP(n,1,5,1) = n21_av(:, 19);
AMP(n,2,5,1) = n22_av(:, 19);
AMP(n,3,5,1) = n23_av(:, 19);
AMP(n,4,5,1) = n24_av(:, 19);
AMP(n,5,5,1) = n25_av(:, 19);
AMP(n,6,5,1) = n26_av(:, 19);
AMP(n,1,5,2) = n21_av(:, 53);
AMP(n,2,5,2) = n22_av(:, 53);
AMP(n,3,5,2) = n23_av(:, 53);
AMP(n,4,5,2) = n24_av(:, 53);
AMP(n,5,5,2) = n25_av(:, 53);
AMP(n,6,5,2) = n26_av(:, 53);
AMP(n,1,5,3) = n21_av(:, 13);
AMP(n,2,5,3) = n22_av(:, 13);
AMP(n,3,5,3) = n23_av(:, 13);
AMP(n,4,5,3) = n24_av(:, 13);
AMP(n,5,5,3) = n25_av(:, 13);
AMP(n,6,5,3) = n26_av(:, 13);
AMP(n,1,5,4) = n21_av(:, 47);
AMP(n,2,5,4) = n22_av(:, 47);
AMP(n,3,5,4) = n23_av(:, 47);
AMP(n,4,5,4) = n24_av(:, 47);
AMP(n,5,5,4) = n25_av(:, 47);
AMP(n,6,5,4) = n26_av(:, 47);
AMP(n,1,5,5) = mean([n21_av(:, 19),n21_av(:, 53),n21_av(:, 13),n21_av(:, 47)],2);
AMP(n,2,5,5) = mean([n22_av(:, 19),n22_av(:, 53),n22_av(:, 13),n22_av(:, 47)],2);
AMP(n,3,5,5) = mean([n23_av(:, 19),n23_av(:, 53),n23_av(:, 13),n23_av(:, 47)],2);
AMP(n,4,5,5) = mean([n24_av(:, 19),n24_av(:, 53),n24_av(:, 13),n24_av(:, 47)],2);
AMP(n,5,5,5) = mean([n25_av(:, 19),n25_av(:, 53),n25_av(:, 13),n25_av(:, 47)],2);
AMP(n,6,5,5) = mean([n26_av(:, 19),n26_av(:, 53),n26_av(:, 13),n26_av(:, 47)],2);

end
%AMP2 = reshape(AMP,Nsub,270); % reshapes: column order = cond(1->6) for each time (1->4), for each electrode (1->5).
AMP2 = reshape(AMP(:,:,:,:),Nsub,150);
AMPLEP = reshape(AMP(:,:,4:5,1),Nsub,12);
%save amplitudes3.mat AMP2
save amplitudes4.mat AMP2
save amplitudesLEP4.mat AMPLEP
P2LAT = squeeze(LATEN(:,:,1));
N2LAT = squeeze(LATEN(:,:,2));
save P2latencies4.mat P2LAT
save N2latencies4.mat N2LAT