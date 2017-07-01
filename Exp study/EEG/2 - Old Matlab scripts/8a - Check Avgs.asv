%Check Averages
clear all

elec_avg_early = [20 51 54];
elec_avg_mid = [15 20 49];
elec_avg_late = [49 50 48 39 14 15 54 51 20]%[14 15 49];
elec_avg_n2a = [14 15 49];
elec_avg_n2b = [13 19 47];
elec_avg_p2 = [15 20 51];

t=-4000:2:1498;

load('4_HE_grand_avg.mat');
grand_avg = iirfilt(grand_avg,500,0,20);
grand_avg = blcorrect4(grand_avg, 250);
plot(t, -grand_avg(15,:), 'b:')
hold on

load('5_HE_grand_avg.mat');
grand_avg = iirfilt(grand_avg,500,0,20);
grand_avg = blcorrect4(grand_avg, 250);
plot(t, -grand_avg(15,:), 'g:')
hold on

load('6_HE_grand_avg.mat');
grand_avg = iirfilt(grand_avg,500,0,20);
grand_avg = blcorrect4(grand_avg, 250);
plot (t, -mean(grand_avg(elec_avg_late,:),1), 'k','LineWidth',2)
hold on

c1 = load('1_HE_grand_avg.mat');
c1.grand_avg = iirfilt(c1.grand_avg,500,0,20);
c2 = load('2_HE_grand_avg.mat');
c2.grand_avg = iirfilt(c2.grand_avg,500,0,20);
c3 = load('3_HE_grand_avg.mat');
c3.grand_avg = iirfilt(c3.grand_avg,500,0,20);
UNC = (c1.grand_avg + c2.grand_avg + c3.grand_avg) / 3;
%plot (t, -UNC(51,:), 'm:')
hold on

c4 = load('4_HE_grand_avg.mat');
c4.grand_avg = iirfilt(c4.grand_avg,500,0,20);
c5 = load('5_HE_grand_avg.mat');
c5.grand_avg = iirfilt(c5.grand_avg,500,0,20);
c6 = load('6_HE_grand_avg.mat');
c6.grand_avg = iirfilt(c6.grand_avg,500,0,20);
CER = (c4.grand_avg + c5.grand_avg + c6.grand_avg) / 3;
SUBT = (c6.grand_avg - c4.grand_avg);
%plot (t, -CER(51,:), 'k:')
hold on

TOT = (CER + UNC) / 2;
HV_TOT = TOT;
HV_high = (c3.grand_avg + c6.grand_avg) / 2;
plot (t, -mean(TOT(elec_avg_late,:),1), 'k','LineWidth',2)
hold on

c6 = load('6_HE_grand_avg.mat');
c6.grand_avg = iirfilt(c6.grand_avg,500,0,20);
[ind peaks] = findpeaks(c6.grand_avg(15, 2100:2400));
P2ind = 2100+ind(find(peaks==max(peaks)));
N2ind = 2100+ind(find(peaks==min(peaks)));
figure
topoplot(mean(c6.grand_avg([1:2 4:30 33:64], 1750:2000),2), 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-5 5]), colorbar
figure
topoplot(mean(c6.grand_avg([1:2 4:30 33:64], 750:1000),2), 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-5 5]), colorbar
%topoplot(mean(TOT([1:2 4:30 33:64], 1250:1500),2), 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-5 5]), colorbar
figure
topoplot(mean(c6.grand_avg([1:2 4:30 33:64], P2ind-5:P2ind+5),2), 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-5 5]), colorbar
%topoplot(mean(TOT([1:2 4:30 33:64], N2ind-5:N2ind+5),2), 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-5 5]), colorbar


%check latency
Y = zeros(1,750);
Y(559:561) = 5;
plot (t, -Y(:,:), 'm')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

%t=-2000:749;

load('4_FM_grand_avg.mat');
plot(t, -grand_avg(49,:), 'b')
hold on

load('5_FM_grand_avg.mat');
plot(t, -grand_avg(49,:), 'g')
hold on

load('6_FM_grand_avg.mat');
grand_avg = iirfilt(grand_avg,500,0,20);
grand_avg = blcorrect4(grand_avg, 250);
plot (t, -mean(grand_avg(elec_avg_late,:),1), 'k--','LineWidth',2)
hold on

c1 = load('1_FM_grand_avg.mat');
c1.grand_avg = iirfilt(c1.grand_avg,500,0,20);
c2 = load('2_FM_grand_avg.mat');
c2.grand_avg = iirfilt(c2.grand_avg,500,0,20);
c3 = load('3_FM_grand_avg.mat');
c3.grand_avg = iirfilt(c3.grand_avg,500,0,20);
UNC = (c1.grand_avg + c2.grand_avg + c3.grand_avg) / 3;
%plot (t, -UNC(51,:), 'm')
hold on

c4 = load('4_FM_grand_avg.mat');
c4.grand_avg = iirfilt(c4.grand_avg,500,0,20);
c5 = load('5_FM_grand_avg.mat');
c5.grand_avg = iirfilt(c5.grand_avg,500,0,20);
c6 = load('6_FM_grand_avg.mat');
CER = (c4.grand_avg + c5.grand_avg + c6.grand_avg) / 3;
SUBT = (c6.grand_avg - c4.grand_avg);
%plot (t, -CER(51,:), 'k')
%hold on

TOT = (CER + UNC) / 2;
FM_TOT = TOT;
FM_high = (c3.grand_avg + c6.grand_avg) / 2;
plot (t, -mean(TOT(elec_avg_late,:),1), 'k--','LineWidth',2)
hold on

c6 = load('6_FM_grand_avg.mat');
c6.grand_avg = iirfilt(c6.grand_avg,500,0,20);
[ind peaks] = findpeaks(c6.grand_avg(15, 2100:2400));
P2ind = 2100+ind(find(peaks==max(peaks)));
N2ind = 2100+ind(find(peaks==min(peaks)));
figure
topoplot(mean(c6.grand_avg([1:2 4:30 33:64], 1750:2000),2), 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-5 5]), colorbar
figure
topoplot(mean(c6.grand_avg([1:2 4:30 33:64], 750:1000),2), 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-5 5]), colorbar
%topoplot(mean(TOT([1:2 4:30 33:64], 1250:1500),2), 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-5 5]), colorbar
figure
topoplot(mean(c6.grand_avg([1:2 4:30 33:64], P2ind-5:P2ind+5),2), 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-5 5]), colorbar
%topoplot(mean(TOT([1:2 4:30 33:64], N2ind-5:N2ind+5),2), 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-5 5]), colorbar



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

%t=-2000:749;

load('4_OA_grand_avg.mat');
plot(t, -grand_avg(49,:), 'b--')
hold on

load('5_OA_grand_avg.mat');
plot(t, -grand_avg(49,:), 'g--')
hold on

load('6_OA_grand_avg.mat');
grand_avg = iirfilt(grand_avg,500,0,20);
grand_avg = blcorrect4(grand_avg, 250);
plot (t, -mean(grand_avg(elec_avg_late,:),1), 'k:','LineWidth',2)
hold on

c1 = load('1_OA_grand_avg.mat');
c1.grand_avg = iirfilt(c1.grand_avg,500,0,20);
c2 = load('2_OA_grand_avg.mat');
c2.grand_avg = iirfilt(c2.grand_avg,500,0,20);
c3 = load('3_OA_grand_avg.mat');
UNC = (c1.grand_avg + c2.grand_avg + c3.grand_avg) / 3;
%plot (t, -UNC(51,:), 'm')
hold on

c4 = load('4_OA_grand_avg.mat');
c4.grand_avg = iirfilt(c4.grand_avg,500,0,20);
c5 = load('5_OA_grand_avg.mat');
c5.grand_avg = iirfilt(c5.grand_avg,500,0,20);
c6 = load('6_OA_grand_avg.mat');
CER = (c4.grand_avg + c5.grand_avg + c6.grand_avg) / 3;
SUBT = (c6.grand_avg - c4.grand_avg);
%plot (t, -CER(51,:), 'k')
%hold on
TOT = (CER + UNC) / 2;
OA_TOT = TOT;
OA_high = (c3.grand_avg + c6.grand_avg) / 2;
plot (t, -mean(TOT(elec_avg_late,:),1), 'k:','LineWidth',2)
hold on

c6 = load('6_OA_grand_avg.mat');
c6.grand_avg = iirfilt(c6.grand_avg,500,0,20);
[ind peaks] = findpeaks(c6.grand_avg(15, 2100:2400));
P2ind = 2100+ind(find(peaks==max(peaks)));
N2ind = 2100+ind(find(peaks==min(peaks)));
figure
topoplot(mean(c6.grand_avg([1:2 4:30 33:64], 1750:2000),2), 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-5 5]), colorbar
figure
topoplot(mean(c6.grand_avg([1:2 4:30 33:64], 750:1000),2), 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-5 5]), colorbar
%topoplot(mean(TOT([1:2 4:30 33:64], 1250:1500),2), 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-5 5]), colorbar
figure
topoplot(mean(c6.grand_avg([1:2 4:30 33:64], P2ind-5:P2ind+5),2), 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-5 5]), colorbar
%topoplot(mean(TOT([1:2 4:30 33:64], N2ind-5:N2ind+5),2), 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-5 5]), colorbar


%GRAND_TOT = (HV_TOT(:,:) + FM_TOT(:,:) + OA_TOT(:,:))/3;
GRAND_TOT = (HV_high(:,:) + FM_high(:,:) + OA_high(:,:))/3;
[ind peaks] = findpeaks(GRAND_TOT(15, 2100:2400));
P2ind = 2100+ind(find(peaks==max(peaks)));
N2ind = 2100+ind(find(peaks==min(peaks)));
topoplot(mean(GRAND_TOT([1:2 4:30 33:64], 1750:2000),2), 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-5 5]), colorbar
topoplot(mean(GRAND_TOT([1:2 4:30 33:64], 750:1000),2), 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-5 5]), colorbar
topoplot(mean(GRAND_TOT([1:2 4:30 33:64], 250:500),2), 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-5 5]), colorbar
topoplot(mean(GRAND_TOT([1:2 4:30 33:64], P2ind-5:P2ind+5),2), 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-5 5]), colorbar
topoplot(mean(GRAND_TOT([1:2 4:30 33:64], N2ind-5:N2ind+5),2), 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-5 5]), colorbar

%late SPN
%load('6_OA_grand_avgP2.mat');
topoplot(mean(TOT([1:2 4:30 33:64], 432:500),2), 'C:\Documents and Settings\All Users\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-5 5]), colorbar

%base
topoplot(mean(TOT([1:2 4:30 33:64], 62:125),2), 'C:\Documents and Settings\All Users\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-5 5]), colorbar

%early
topoplot(mean(TOT([1:2 4:30 33:64], 188:250),2), 'C:\Documents and Settings\All Users\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-5 5]), colorbar

%mid
topoplot(mean(TOT([1:2 4:30 33:64], 313:375),2), 'C:\Documents and Settings\All Users\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-5 5]), colorbar

%late
topoplot(mean(TOT([1:2 4:30 33:64], 1750:2000),2), 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-5 5]), colorbar

%P2 HE
topoplot(mean(grand_avg([1:2 4:30 33:64], 559:561),2), 'C:\Documents and Settings\All Users\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-6 6]), colorbar

%P2 FM
topoplot(mean(grand_avg([1:2 4:30 33:64], 561:563),2), 'C:\Documents and Settings\All Users\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-6 6]), colorbar

%N2 HE
topoplot(mean(grand_avg([1:2 4:30 33:64], 542:544),2), 'C:\Documents and Settings\All Users\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-6 6]), colorbar

%N2 FM
topoplot(mean(grand_avg([1:2 4:30 33:64], 539:540),2), 'C:\Documents and Settings\All Users\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-6 6]), colorbar

%N2 OA
topoplot(mean(grand_avg([1:2 4:30 33:64], 548:549),2), 'C:\Documents and Settings\All Users\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-6 6]), colorbar

%VEP 1
topoplot(mean(TOT([1:2 4:30 33:64], 134),2), 'C:\Documents and Settings\All Users\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-6 6]), colorbar

