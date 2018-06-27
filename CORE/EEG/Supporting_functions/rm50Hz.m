function clean = rm50Hz(data,ord,sr,ncomp) % data, order, no components to remove
%ord = order of the data with respect to samples x 

% This example uses DSS to remove power-line noise from EEG data.
% Finds components that are most dominated by 50 Hz and 
% harmonics, regresses them out to obtain clean data.
%
% Uses nt_bias_fft(), nt_dss0(), nt_tsr().

% use DSS to isolate 50 Hz & harmonics
disp('50 Hz & harmonics DSS...');

% re-order data
data = permute(data,ord);

% covariance matrices of full band (c0) and filtered to 50 Hz & harmonics (c1)
[c0,c1]=nt_bias_fft(data,[50, 100, 150]/sr, 512);

% DSS matrix
[todss,pwr0,pwr1]=nt_dss0(c0,c1); 
p1=pwr1./pwr0; % score, proportional to power ratio of 50Hz & harmonics to full band

% DSS components
z=nt_mmat(data,todss);

% first components are most dominated by 50Hz & harmonics
NREMOVE=ncomp;
clean=nt_tsr(data,z(:,1:NREMOVE,:)); % regress them out

% plot bias score
figure(1); clf; set(gcf,'color', [1 1 1]);
plot(p1, '.-'); xlabel('component'); ylabel('score'); title('DSS score');

% plot spectra of DSS components
figure(2); clf; set(gcf,'color', [1 1 1]);
nt_spect_plot2(nt_normcol(z(:,1:30,:)),512,[],[],sr);
title('spectra of first 30 DSS components'); ylabel('component')

% plot spectra of data before and after removal of power line components
figure(3); clf; set(gcf,'color', [1 1 1]);
subplot 121;
nt_spect_plot(data,1024,[],[],sr); 
set(gca,'ygrid','on');
hold on
nt_spect_plot(clean,1024,[],[],sr);
nt_colorlines([],[3 1]);
title('power spectra, average over channels');
legend('before','after'); legend boxoff
set(gca,'ygrid','on');
subplot 122
nt_spect_plot(data-clean,1024,[],[],sr);
title('noise power (removed)');
set(gca,'ygrid','on');


clean = ipermute(clean,ord);


