%extract data from Cz, FCz, average of (PO7, P5, P7), average of (PO8, P6, P8)

clear all

subjects = {'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13' 'H14', 'H15', 'H16', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16'};
   
Nsub = length(subjects);

AMP = zeros(Nsub,6,2,4); % subject, condition, time bins, electrode
LATEN = zeros(Nsub,6,2);

fnames={'_1_avg_ca.mat';
    '_2_avg_ca.mat';
    '_3_avg_ca.mat';
    '_4_avg_ca.mat';
    '_5_avg_ca.mat';
    '_6_avg_ca.mat';
  };

for n = 1:Nsub
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

[XMAX,IMAX,XMIN,IMIN] = extrema(Cz(:,120:160));%%search for extrema on internet - peak detection prog, xmax ampli imax latency
LATp = IMAX(1)+119; %% change to be the same as the 1st time point to zero everything
LATn = IMIN(1)+119;

[XMAX,IMAX,XMIN,IMIN] = extrema(mdata1.avg(15,LATp-8:LATp+8));%%each cond seperatly +/-64 ms 15 is CZ
LATp1 = IMAX(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata1.avg(15,LATn-8:LATn+8));
LATn1 = IMIN(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata2.avg(15,LATp-8:LATp+8));
LATp2 = IMAX(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata2.avg(15,LATn-8:LATn+8));
LATn2 = IMIN(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata3.avg(15,LATp-8:LATp+8));
LATp3 = IMAX(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata3.avg(15,LATn-8:LATn+8));
LATn3 = IMIN(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata4.avg(15,LATp-8:LATp+8));
LATp4 = IMAX(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata4.avg(15,LATn-8:LATn+8));
LATn4 = IMIN(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata5.avg(15,LATp-8:LATp+8));
LATp5 = IMAX(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata5.avg(15,LATn-8:LATn+8));
LATn5 = IMIN(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata6.avg(15,LATp-8:LATp+8));
LATp6 = IMAX(1);
[XMAX,IMAX,XMIN,IMIN] = extrema(mdata6.avg(15,LATn-8:LATn+8));
LATn6 = IMIN(1);

%Extract and save data

p1=mdata1.avg(:,(LATp-9+LATp1)-1:(LATp-9+LATp1)+1)';%P for condtion 1 
p2=mdata2.avg(:,(LATp-9+LATp2)-1:(LATp-9+LATp2)+1)';%P for condition 2
p3=mdata3.avg(:,(LATp-9+LATp3)-1:(LATp-9+LATp3)+1)';
p4=mdata4.avg(:,(LATp-9+LATp4)-1:(LATp-9+LATp4)+1)';
p5=mdata5.avg(:,(LATp-9+LATp5)-1:(LATp-9+LATp5)+1)';
p6=mdata6.avg(:,(LATp-9+LATp6)-1:(LATp-9+LATp6)+1)';
p1_av=mean(p1);
p2_av=mean(p2);
p3_av=mean(p3);
p4_av=mean(p4);
p5_av=mean(p5);
p6_av=mean(p6);
%c_p2_av=[p24_av,p25_av,p26_av];
%c_p2_av=mean(c_p2_av,2);
%u_p2_av=[p21_av,p22_av,p23_av];
%u_p2_av=mean(u_p2_av,2);

eval(['save ' subject '_1_vp_av.dat'  ' p1_av /ASCII'])%all electrodes in 1 file - source analysis
eval(['save ' subject '_2_vp_av.dat'  ' p2_av /ASCII'])
eval(['save ' subject '_3_vp_av.dat'  ' p3_av /ASCII'])
eval(['save ' subject '_4_vp_av.dat'  ' p4_av /ASCII'])
eval(['save ' subject '_5_vp_av.dat'  ' p5_av /ASCII'])
eval(['save ' subject '_6_vp_av.dat'  ' p6_av /ASCII'])

LATEN(n,1,1) = ((LATp-9+LATp1)*8)-4000;
LATEN(n,2,1) = ((LATp-9+LATp2)*8)-4000;
LATEN(n,3,1) = ((LATp-9+LATp3)*8)-4000;
LATEN(n,4,1) = ((LATp-9+LATp4)*8)-4000;
LATEN(n,5,1) = ((LATp-9+LATp5)*8)-4000;
LATEN(n,6,1) = ((LATp-9+LATp6)*8)-4000;

AMP(n,1,1,1) = p1_av(:, 15);
AMP(n,2,1,1) = p2_av(:, 15);
AMP(n,3,1,1) = p3_av(:, 15);
AMP(n,4,1,1) = p4_av(:, 15);
AMP(n,5,1,1) = p5_av(:, 15);
AMP(n,6,1,1) = p6_av(:, 15);
AMP(n,1,1,2) = p1_av(:, 20);
AMP(n,2,1,2) = p2_av(:, 20);
AMP(n,3,1,2) = p3_av(:, 20);
AMP(n,4,1,2) = p4_av(:, 20);
AMP(n,5,1,2) = p5_av(:, 20);
AMP(n,6,1,2) = p6_av(:, 20);
AMP(n,1,1,3) = mean([p1_av(:,11),p1_av(:,45),p1_av(:,46)],2);
AMP(n,2,1,3) = mean([p2_av(:,11),p2_av(:,45),p2_av(:,46)],2);
AMP(n,3,1,3) = mean([p3_av(:,11),p3_av(:,45),p3_av(:,46)],2);
AMP(n,4,1,3) = mean([p4_av(:,11),p4_av(:,45),p4_av(:,46)],2);
AMP(n,5,1,3) = mean([p5_av(:,11),p5_av(:,45),p5_av(:,46)],2);
AMP(n,6,1,3) = mean([p6_av(:,11),p6_av(:,45),p6_av(:,46)],2);
AMP(n,1,1,4) = mean([p1_av(:,7),p1_av(:,40),p1_av(:,44)],2);
AMP(n,2,1,4) = mean([p2_av(:,7),p2_av(:,40),p2_av(:,44)],2);
AMP(n,3,1,4) = mean([p3_av(:,7),p3_av(:,40),p3_av(:,44)],2);
AMP(n,4,1,4) = mean([p4_av(:,7),p4_av(:,40),p4_av(:,44)],2);
AMP(n,5,1,4) = mean([p5_av(:,7),p5_av(:,40),p5_av(:,44)],2);
AMP(n,6,1,4) = mean([p6_av(:,7),p6_av(:,40),p6_av(:,44)],2);

n1=mdata1.avg(:,(LATn-9+LATn1)-1:(LATn-9+LATn1)+1)';%n2 for condtion 1 
n2=mdata2.avg(:,(LATn-9+LATn2)-1:(LATn-9+LATn2)+1)';%n2 for condition 2
n3=mdata3.avg(:,(LATn-9+LATn3)-1:(LATn-9+LATn3)+1)';
n4=mdata4.avg(:,(LATn-9+LATn4)-1:(LATn-9+LATn4)+1)';
n5=mdata5.avg(:,(LATn-9+LATn5)-1:(LATn-9+LATn5)+1)';
n6=mdata6.avg(:,(LATn-9+LATn6)-1:(LATn-9+LATn6)+1)';
n1_av=mean(n1);
n2_av=mean(n2);
n3_av=mean(n3);
n4_av=mean(n4);
n5_av=mean(n5);
n6_av=mean(n6);
%c_n2_av=[n24_av,n25_av,n26_av];
%c_n2_av=mean(c_n2_av,2);
%u_n2_av=[n21_av,n22_av,n23_av];
%u_n2_av=mean(u_n2_av,2);

eval(['save ' subject '_1_vn_av.dat'  ' n1_av /ASCII'])%all electrodes in 1 file - source analysis
eval(['save ' subject '_2_vn_av.dat'  ' n2_av /ASCII'])
eval(['save ' subject '_3_vn_av.dat'  ' n3_av /ASCII'])
eval(['save ' subject '_4_vn_av.dat'  ' n4_av /ASCII'])
eval(['save ' subject '_5_vn_av.dat'  ' n5_av /ASCII'])
eval(['save ' subject '_6_vn_av.dat'  ' n6_av /ASCII'])

LATEN(n,1,2) = ((LATn-9+LATn1)*8)-4000;
LATEN(n,2,2) = ((LATn-9+LATn2)*8)-4000;
LATEN(n,3,2) = ((LATn-9+LATn3)*8)-4000;
LATEN(n,4,2) = ((LATn-9+LATn4)*8)-4000;
LATEN(n,5,2) = ((LATn-9+LATn5)*8)-4000;
LATEN(n,6,2) = ((LATn-9+LATn6)*8)-4000;

AMP(n,1,2,1) = n1_av(:, 15);
AMP(n,2,2,1) = n2_av(:, 15);
AMP(n,3,2,1) = n3_av(:, 15);
AMP(n,4,2,1) = n4_av(:, 15);
AMP(n,5,2,1) = n5_av(:, 15);
AMP(n,6,2,1) = n6_av(:, 15);
AMP(n,1,2,2) = n1_av(:, 20);
AMP(n,2,2,2) = n2_av(:, 20);
AMP(n,3,2,2) = n3_av(:, 20);
AMP(n,4,2,2) = n4_av(:, 20);
AMP(n,5,2,2) = n5_av(:, 20);
AMP(n,6,2,2) = n6_av(:, 20);
AMP(n,1,2,3) = mean([n1_av(:,11),n1_av(:,45),n1_av(:,46)],2);
AMP(n,2,2,3) = mean([n2_av(:,11),n2_av(:,45),n2_av(:,46)],2);
AMP(n,3,2,3) = mean([n3_av(:,11),n3_av(:,45),n3_av(:,46)],2);
AMP(n,4,2,3) = mean([n4_av(:,11),n4_av(:,45),n4_av(:,46)],2);
AMP(n,5,2,3) = mean([n5_av(:,11),n5_av(:,45),n5_av(:,46)],2);
AMP(n,6,2,3) = mean([n6_av(:,11),n6_av(:,45),n6_av(:,46)],2);
AMP(n,1,2,4) = mean([n1_av(:,7),n1_av(:,40),n1_av(:,44)],2);
AMP(n,2,2,4) = mean([n2_av(:,7),n2_av(:,40),n2_av(:,44)],2);
AMP(n,3,2,4) = mean([n3_av(:,7),n3_av(:,40),n3_av(:,44)],2);
AMP(n,4,2,4) = mean([n4_av(:,7),n4_av(:,40),n4_av(:,44)],2);
AMP(n,5,2,4) = mean([n5_av(:,7),n5_av(:,40),n5_av(:,44)],2);
AMP(n,6,2,4) = mean([n6_av(:,7),n6_av(:,40),n6_av(:,44)],2);

end
AMP2 = reshape(AMP,Nsub,48); % reshapes: column order = cond(1->6) for each time (1->4), for each electrode (1->5).
save amplitudesVEP.mat AMP2
PLAT = squeeze(LATEN(:,:,1));
NLAT = squeeze(LATEN(:,:,2));
save vPlatencies.mat PLAT
save vNlatencies.mat NLAT