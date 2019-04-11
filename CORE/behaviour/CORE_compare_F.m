
% load('CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061503_HC_it18.mat')
% Fh=F
% load('CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_CPRS_it18.mat')
% Fc=F
clear all
cd('C:\Data\CORE\behaviour\hgf\fitted')
load('D_fit_r1_it10_n30r.mat')
Fh=GS(2).F;
Fc=GS(1).F;
    
figure
plot(Fc,'r'); hold on; plot(Fh,'b')
thresh = abs((Fc(1)+Fh(1))/2 * 0.01);
for i = 2:length(Fh)
    Fdiff(1,i) = abs(Fh(i)-Fh(i-1));
    Fdiff(2,i) = abs(Fc(i)-Fc(i-1));
end
figure
plot(Fdiff(2,:),'r'); hold on; plot(Fdiff(1,:),'b')
hold on
line([0 length(Fh)],[thresh thresh])