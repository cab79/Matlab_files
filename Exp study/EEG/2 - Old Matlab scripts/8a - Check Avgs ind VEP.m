%Check Averages
clear all

subjects = {'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13' 'H14', 'H15', 'H16', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'OA1', 'OA2', 'OA4', 'OA5', 'OA6', 'OA7', 'OA8', 'OA9', 'OA10', 'OA11', 'OA12'};

nSub = 1;
subject = subjects(nSub);
subject = char(subject);

t=-500:250;

%load([subject '_1_avg.mat']);
%plot(t, -avg(15,:), 'b:');
%hold on

%load([subject '_2_avg.mat']);
%plot(t, -avg(15,:), 'g:');
%hold on

%load([subject '_3_avg.mat']);
%plot(t, -avg(15,:), 'r:');
%hold on

load([subject '_4_avg_ca.mat']);
plot(t, -avg(15,:), 'k');
hold on

load([subject '_5_avg_ca.mat']);
plot(t, -avg(15,:), 'k');
hold on

load([subject '_6_avg_ca.mat']);
plot(t, -avg(15,:), 'k');
hold on

%check latency
load vNlatencies
load vPlatencies
%load amplitudesLEP
load amplitudesVEP
N = (NLAT+4000)/8;
P = (PLAT+4000)/8;
Y = zeros(1,length(t));
Z = zeros(1,length(t));

%Y(P2(nSub,1)) = AMPLEP(nSub,1);
%Y(P2(nSub,2)) = AMPLEP(nSub,2);
%Y(P2(nSub,3)) = AMPLEP(nSub,3);
Y(P(nSub,4)) = AMP2(nSub,4);
plot (t, -Y(:,:), 'y.')
hold on
Y(P(nSub,5)) = AMP2(nSub,5);
plot (t, -Y(:,:), 'y.')
hold on
Y(P(nSub,6)) = AMP2(nSub,6);
plot (t, -Y(:,:), 'y.')
hold on

%Z(N2(nSub,1)) = AMPLEP(nSub,7);
%Z(N2(nSub,2)) = AMPLEP(nSub,8);
%Z(N2(nSub,3)) = AMPLEP(nSub,9);
Z(N(nSub,4)) = AMP2(nSub,10);
plot (t, -Z(:,:), 'm.')
hold on
Z(N(nSub,5)) = AMP2(nSub,11);
plot (t, -Z(:,:), 'm.')
hold on
Z(N(nSub,6)) = AMP2(nSub,12);
plot (t, -Z(:,:), 'm.')
hold on