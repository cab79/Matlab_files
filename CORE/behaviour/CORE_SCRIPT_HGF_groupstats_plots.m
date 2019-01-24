% extract, tabulate and save parameters and summary stats of trajectories

close all
clear all
dbstop if error
restoredefaultpath

% add toolbox paths
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')
addpath('C:\Data\Matlab\raincloud_plots');
addpath('C:\Data\Matlab\Violinplot-Matlab-master');
addpath('C:\Data\Matlab\cbrewer'); % (https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab)

% file info
S.path.hgf = 'C:\Data\CORE\behaviour\hgf'; 

fnames = {
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20180821T134505.mat'; % FITTED %alpha prior = 0.5, cp50 effects
% %'CORE_fittedparameters_percmodel12_respmodel20_fractrain0_20181021T093241.mat';
% %'CORE_fittedparameters_percmodel12_bayesopt_20181019T083824.mat'; % BAYESOPT
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181023T212308.mat'; %alpha prior = 1
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181023T202826.mat'; %alpha prior = 0.2, cp50 effects
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181024T222027.mat'; %fitted with bayesopt priors
% 
% % comparison of prior variances: bayesopt omegas
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181026T171952.mat'; %fitted with bayesopt priors, var025 
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181026T171959.mat'; %fitted with bayesopt priors, var05
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181026T171938.mat'; %fitted with bayesopt priors, var1, cp50 effects
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181026T172007.mat'; %fitted with bayesopt priors, var2, cp50 effects
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181026T172024.mat'; %fitted with bayesopt priors, var4
% 
% % prior on expected uncertainty = -2
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181028T082701.mat'; %fitted with bayesopt priors, var05
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181028T082725.mat'; %fitted with bayesopt priors, var1
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181028T082749.mat'; %fitted with bayesopt priors, var2
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181028T082828.mat'; %fitted with bayesopt priors, var4
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181028T082859.mat'; %fitted with bayesopt priors, var8
% 
% % alpha prior = 0.01
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181101T211928.mat'; %fitted with bayesopt priors, var05
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181101T211935.mat'; %fitted with bayesopt priors, var1
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181101T211950.mat'; %fitted with bayesopt priors, var2
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181101T211959.mat'; %fitted with bayesopt priors, var4
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181101T212006.mat'; %fitted with bayesopt priors, var8 Significant Dau MVA
% % alpha prior = 1
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181102T192113.mat'; %fitted with bayesopt priors, var05 Dau MVA 0.34 cond
%  'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181102T192126.mat'; %fitted with bayesopt priors, var1 Significant Dau MVA **
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181102T192137.mat'; %fitted with bayesopt priors, var2 Dau MVA 0.31 cond
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181102T192151.mat'; %fitted with bayesopt priors, var4 Dau MVA 0.34 cond
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181102T192201.mat'; %fitted with bayesopt priors, var8
% % alpha prior = 0.5
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181104T172839.mat'; %fitted with bayesopt priors, var05
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181104T172830.mat'; %fitted with bayesopt priors, var1
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181104T172850.mat'; %fitted with bayesopt priors, var2 Significant Dau MVA
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181104T172900.mat'; %fitted with bayesopt priors, var4
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181104T172909.mat'; %fitted with bayesopt priors, var8
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181106T212138.mat'; %fitted with bayesopt priors, var10 Dau MVA 0.27 cond, oddball cond
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181106T211643.mat'; %fitted with bayesopt priors, var100
% 
% % ommu2 = -2
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181107T214827.mat'; % ommu3=-3 Significant Dau MVA, oddball cond
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181107T214552.mat'; % ommu3=-4
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181107T214842.mat'; % ommu3=-5
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181107T214441.mat'; % ommu3=-6
% % ommu2 = -6
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181110T095244.mat'; % ommu3=-3
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181109T200448.mat'; % ommu3=-4
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181109T200514.mat'; % ommu3=-5
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181109T200546.mat'; % ommu3=-6
% 
% % alphas with minimum omega priors
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181114T173639.mat'; % alpha 0.1
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181114T173746.mat'; % alpha 0.2
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181114T173907.mat'; % alpha 0.4 % subtracted cond effects
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181114T173922.mat'; % alpha 0.8
% 
% %alphas with minimum omega -5
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181116T115736.mat'; % alpha 0.1 Dau MVA 0.34 cond
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181116T115701.mat'; % alpha 0.2
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181116T115619.mat'; % alpha 0.4
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181116T115537.mat'; % alpha 0.8
% 
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181119T111012.mat' % 'GBM_config_alpha01_var4'
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181119T111020.mat' % 'GBM_config_alpha02_var4'
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181119T111029.mat' % 'GBM_config_alpha04_var4'
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181119T111044.mat' % 'GBM_config_alpha08_var4'
% 
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181117T214940.mat' % 'GBM_config_alpha01_var2'
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181117T214948.mat' % 'GBM_config_alpha02_var2' , cp50 effects
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181117T214959.mat' % 'GBM_config_alpha04_var2' , cp50 effects
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181117T215010.mat' % 'GBM_config_alpha08_var2'
% 
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181119T211705.mat' % 'GBM_config_alpha01_var1'
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181119T211717.mat' % 'GBM_config_alpha02_var1' , cp50 effects
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181119T211727.mat' % 'GBM_config_alpha04_var1'
%  'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181119T211738.mat' % 'GBM_config_alpha08_var1' % subtracted cond effects
% 
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181121T004327.mat' % 'GBM_config_alpha01_var05'
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181120T235241.mat' % 'GBM_config_alpha02_var05'
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181121T024806.mat' % 'GBM_config_alpha04_var05' Significant Dau MVA
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181121T090841.mat' % 'GBM_config_alpha08_var05'
% 
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181121T151444.mat' % 'GBM_config_alpha01_var4_bo'
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181121T232554.mat' % 'GBM_config_alpha02_var4_bo'
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181121T223551.mat' % 'GBM_config_alpha04_var4_bo'
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181121T124130.mat' % 'GBM_config_alpha08_var4_bo'
% 
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181122T211624.mat' % 'GBM_config_alpha01_var2_bo'
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181122T211634.mat' % 'GBM_config_alpha02_var2_bo' oddball cond
%  'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181122T211720.mat' % 'GBM_config_alpha04_var2_bo'
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181122T213340.mat' % 'GBM_config_alpha08_var2_bo'
% 
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181123T111508.mat' % 'GBM_config_alpha01_var1_bo'
%  'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181123T111519.mat' % 'GBM_config_alpha02_var1_bo' Dau MVA 0.27 cond
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181123T111527.mat' % 'GBM_config_alpha04_var1_bo'
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181123T111534.mat' % 'GBM_config_alpha08_var1_bo' Significant Dau MVA ** Also 0.34 cond
% 
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181124T142822.mat' % 'GBM_config_alpha01_var05_bo'
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181124T142836.mat' % 'GBM_config_alpha02_var05_bo'
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181124T142849.mat' % 'GBM_config_alpha04_var05_bo'
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181124T142907.mat' % 'GBM_config_alpha08_var05_bo' oddball cond, cp50 effects
% 
% % one alpha
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181126T085335.mat' % 'GBM_config_alpha01_var4'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181126T085344.mat' % 'GBM_config_alpha02_var4'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181126T085414.mat' % 'GBM_config_alpha04_var4'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181126T084833.mat' % 'GBM_config_alpha01_var4_bo'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181126T084924.mat' % 'GBM_config_alpha02_var4_bo'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181126T084948.mat' % 'GBM_config_alpha04_var4_bo'
% % 
% % one alpha BO
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181127T094048.mat'; % 'GBM_config_alphaBO_var4_bo'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181127T094128.mat'; % 'GBM_config_alpha2BO_var4_bo'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181127T094159.mat'; % 'GBM_config_alpha4BO_var4_bo'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181127T094246.mat'; % 'GBM_config_alpha8BO_var4_bo'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181128T091331.mat'; % 'GBM_config_alphaBO_var2_bo'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181128T091401.mat'; % 'GBM_config_alpha2BO_var2_bo'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181128T091427.mat'; % 'GBM_config_alpha4BO_var2_bo'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181128T091450.mat'; % 'GBM_config_alpha8BO_var2_bo'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181128T091518.mat'; % 'GBM_config_alphaBO_var1_bo'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181128T091544.mat'; % 'GBM_config_alpha2BO_var1_bo'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181128T091611.mat'; % 'GBM_config_alpha4BO_var1_bo'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181128T091646.mat'; % 'GBM_config_alpha8BO_var1_bo'
% 
% % S.path.hgf = 'C:\Data\CORE\behaviour\hgf\sim';
% % fname = 'CORE_sim_percmodel12_20181019T101036.mat'; % Al*2, om2+2, om3-2
% % fname = 'CORE_sim_percmodel12_20181019T102209.mat'; % Al*2
% % fname = 'CORE_sim_percmodel12_20181019T102418.mat'; % om2+2
% % fname = 'CORE_sim_percmodel12_20181019T102607.mat'; % om3-2
% % fname = 'CORE_sim_percmodel12_20181019T102839.mat'; % om2+2, om3-2
% % fname = 'CORE_sim_percmodel12_20181019T103111.mat'; % Al*2, om2+2
% % fname = 'CORE_sim_percmodel12_20181019T103323.mat'; % Al*2, om3-2
% % 
% % one alpha BO
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181207T112500.mat'; % 'GBM_config_alpha1.5BO_var1.5_bo'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181128T091401.mat'; % 'GBM_config_alpha2BO_var2_bo'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181207T112310.mat'; % 'GBM_config_alpha3BO_var3_bo'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181127T094159.mat'; % 'GBM_config_alpha4BO_var4_bo'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181207T112129.mat'; % 'GBM_config_alpha5BO_var5_bo'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181207T111915.mat'; % 'GBM_config_alpha6BO_var6_bo'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181207T111853.mat'; % 'GBM_config_alpha7BO_var7_bo'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181207T110417.mat'; % 'GBM_config_alpha8BO_var8_bo'
% % 
% % %logvar
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181214T070954.mat'; %'GBM_config_alpha05BO_var4_bo'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181214T071046.mat'; %'GBM_config_alpha1BO_var4_bo'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181214T071111.mat' %'GBM_config_alpha1.1BO_var4_bo'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181214T071147.mat'; %'GBM_config_alpha1.5BO_var4_bo'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181214T071217.mat'; %'GBM_config_alpha2BO_var4_bo' % signifcant Dau using LDA
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181214T071303.mat'; %'GBM_config_alpha4BO_var4_bo'
% % 
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181216T093231.mat'; %'GBM_config_alpha1BO_var4_bo'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181216T093251.mat'; %'GBM_config_alpha1.5BO_var4_bo'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181216T093310.mat'; %'GBM_config_alpha2BO_var4_bo'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181216T093329.mat'; %'GBM_config_alpha4BO_var4_bo'
% % 'CORE_fittedparameters_percmodel3_respmodel2_fractrain0_20181216T093342.mat'; %'GBM_config_alpha8BO_var4_bo'
% % 
% % %resp model 21, one alpha
% % 'CORE_fittedparameters_percmodel3_respmodel21_fractrain0_20181216T093231.mat' %'GBM_config_alpha1BO_var4_bo'
% 
% %4 alpha, 3 resp models
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181219T170032.mat'; %'GBM_config_alpha1BO_var4_bo4'
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181219T212858.mat'; %'GBM_config_alpha2BO_var4_bo4'
% 'CORE_fittedparameters_percmodel12_respmodel21_fractrain0_20181219T170032.mat'; %'GBM_config_alpha1BO_var4_bo4'
% 'CORE_fittedparameters_percmodel12_respmodel21_fractrain0_20181219T212858.mat'; %'GBM_config_alpha2BO_var4_bo4'
% 'CORE_fittedparameters_percmodel12_respmodel12_fractrain0_20181219T170032.mat'; %'GBM_config_alpha1BO_var4_bo4'
% 'CORE_fittedparameters_percmodel12_respmodel12_fractrain0_20181219T212858.mat'; %'GBM_config_alpha2BO_var4_bo4'
% 
% %4 increasing mean on all params, alpha var 1
% 'CORE_fittedparameters_percmodel1205_respmodel2_fractrain0_20181223T102343.mat'; %'GBM_config_05BO_al4'
% 'CORE_fittedparameters_percmodel121_respmodel2_fractrain0_20181223T102343.mat'; %'GBM_config_1BO_al4'
% 'CORE_fittedparameters_percmodel122_respmodel2_fractrain0_20181223T102343.mat'; %'GBM_config_2BO_al4'
% 'CORE_fittedparameters_percmodel124_respmodel2_fractrain0_20181223T102343.mat'; %'GBM_config_4BO_al4'
% 'CORE_fittedparameters_percmodel128_respmodel2_fractrain0_20181223T102343.mat'; %'GBM_config_8BO_al4'
% 'CORE_fittedparameters_percmodel1216_respmodel2_fractrain0_20181223T102343.mat'; %'GBM_config_16BO_al4'
% 'CORE_fittedparameters_percmodel1232_respmodel2_fractrain0_20181223T102343.mat'; %'GBM_config_32BO_al4' 
% 
% %4 increasing mean on all params, alpha var 2
% 'CORE_fittedparameters_percmodel1205_respmodel2_fractrain0_20181224T231959.mat'; %'GBM_config_05BO_al4'
% 'CORE_fittedparameters_percmodel121_respmodel2_fractrain0_20181224T231959.mat'; %'GBM_config_1BO_al4'
% 'CORE_fittedparameters_percmodel122_respmodel2_fractrain0_20181224T231959.mat'; %'GBM_config_2BO_al4'
% 'CORE_fittedparameters_percmodel124_respmodel2_fractrain0_20181224T231959.mat'; %'GBM_config_4BO_al4'
% 'CORE_fittedparameters_percmodel128_respmodel2_fractrain0_20181224T231959.mat'; %'GBM_config_8BO_al4'
% 'CORE_fittedparameters_percmodel1216_respmodel2_fractrain0_20181224T231959.mat'; %'GBM_config_16BO_al4'
% 'CORE_fittedparameters_percmodel1232_respmodel2_fractrain0_20181224T231959.mat'; %'GBM_config_32BO_al4' 

% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181226T153015.mat' % 'GBM_config_1BO_al4_alvar1' max 1, fill 2
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181226T153145.mat' % 'GBM_config_1BO_al4_alvar2' max 1, fill 2
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181226T152946.mat' % 'GBM_config_1BO_al4_alvar1' max 0.2, fill 0
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181226T153159.mat' % 'GBM_config_1BO_al4_alvar2' max 0.2, fill 0
% 'CORE_fittedparameters_percmodel12_respmodel2_fractrain0_20181226T153335.mat' % 'GBM_config_1BO_al4_alvar4' max 1, fill 0

% is model 15 or 16 best?
% 'CORE_fittedparameters_percmodel121_respmodel15_fractrain0_20181227T172306.mat' %'GBM_config_1BO_al4_alvar1'
% 'CORE_fittedparameters_percmodel121_respmodel16_fractrain0_20181227T172306.mat' %'GBM_config_1BO_al4_alvar1'

% % alpha om2/3 increasing
% 'CORE_fittedparameters_percmodel121_respmodel15_fractrain0_20181228T104344.mat'; %'GBM_config_1BO_al4'
% 'CORE_fittedparameters_percmodel122_respmodel15_fractrain0_20181228T104344.mat'; %'GBM_config_2BO_al4'
% 'CORE_fittedparameters_percmodel124_respmodel15_fractrain0_20181228T104344.mat'; %'GBM_config_4BO_al4'
% 'CORE_fittedparameters_percmodel128_respmodel15_fractrain0_20181228T104344.mat'; %'GBM_config_8BO_al4'
% 'CORE_fittedparameters_percmodel1216_respmodel15_fractrain0_20181228T104344.mat'; %'GBM_config_16BO_al4' 
% 'CORE_fittedparameters_percmodel1232_respmodel15_fractrain0_20181228T104344.mat'; %'GBM_config_32BO_al4'
% 'CORE_fittedparameters_percmodel1264_respmodel15_fractrain0_20181228T104344.mat'; %'GBM_config_64BO_al4'
% 
% % alpha increasing
% 'CORE_fittedparameters_percmodel121_respmodel15_fractrain0_20181228T193827.mat'; %'GBM_config_1BO_al4'
% 'CORE_fittedparameters_percmodel122_respmodel15_fractrain0_20181228T193827.mat'; %'GBM_config_2BO_al4'
% 'CORE_fittedparameters_percmodel124_respmodel15_fractrain0_20181228T193827.mat'; %'GBM_config_4BO_al4'
% 'CORE_fittedparameters_percmodel128_respmodel15_fractrain0_20181228T193827.mat'; %'GBM_config_8BO_al4'
% 'CORE_fittedparameters_percmodel1216_respmodel15_fractrain0_20181228T193827.mat'; %'GBM_config_16BO_al4' 
% 'CORE_fittedparameters_percmodel1232_respmodel15_fractrain0_20181228T193827.mat'; %'GBM_config_32BO_al4' 
% 'CORE_fittedparameters_percmodel1264_respmodel15_fractrain0_20181228T193827.mat'; %'GBM_config_64BO_al4' 

% % % alpha variance increasing
% 'CORE_fittedparameters_percmodel121_respmodel15_fractrain0_20181231T111657.mat'; %'GBM_config_BO_al4_alvar1'
% 'CORE_fittedparameters_percmodel122_respmodel15_fractrain0_20181231T111657.mat'; %'GBM_config_BO_al4_alvar2'
% 'CORE_fittedparameters_percmodel124_respmodel15_fractrain0_20181231T111657.mat'; %'GBM_config_BO_al4_alvar4'
% 'CORE_fittedparameters_percmodel128_respmodel15_fractrain0_20181231T111657.mat'; %'GBM_config_BO_al4_alvar8'
% 'CORE_fittedparameters_percmodel1216_respmodel15_fractrain0_20181231T111657.mat'; %'GBM_config_BO_al4_alvar16'  
% 
% % % alpha variance increasing
% 'CORE_fittedparameters_percmodel121_respmodel15_fractrain0_20181231T120153.mat'; %'GBM_config_2BO_al4_alvar1'
% 'CORE_fittedparameters_percmodel122_respmodel15_fractrain0_20181231T120153.mat'; %'GBM_config_2BO_al4_alvar2'
% 'CORE_fittedparameters_percmodel124_respmodel15_fractrain0_20181231T120153.mat'; %'GBM_config_2BO_al4_alvar4'
% 'CORE_fittedparameters_percmodel128_respmodel15_fractrain0_20181231T120153.mat'; %'GBM_config_2BO_al4_alvar8'
% 'CORE_fittedparameters_percmodel1216_respmodel15_fractrain0_20181231T120153.mat'; %'GBM_config_2BO_al4_alvar16'  

% % alpha mean increasing, alvar4
% 'CORE_fittedparameters_percmodel1205_respmodel15_fractrain0_20190101T162309.mat'; %'GBM_config_05BO_al4_alvar4'
% 'CORE_fittedparameters_percmodel1210_respmodel15_fractrain0_20190101T162309.mat'; %'GBM_config_10BO_al4_alvar4'
% 'CORE_fittedparameters_percmodel1215_respmodel15_fractrain0_20190101T162309.mat'; %'GBM_config_15BO_al4_alvar4'
% 'CORE_fittedparameters_percmodel1220_respmodel15_fractrain0_20190101T162309.mat'; %'GBM_config_20BO_al4_alvar4'
% 'CORE_fittedparameters_percmodel1225_respmodel15_fractrain0_20190101T162309.mat'; %'GBM_config_25BO_al4_alvar4' 
% 'CORE_fittedparameters_percmodel1230_respmodel15_fractrain0_20190101T162309.mat'; %'GBM_config_30BO_al4_alvar4'
% 'CORE_fittedparameters_percmodel1235_respmodel15_fractrain0_20190101T162309.mat'; %'GBM_config_35BO_al4_alvar4'
% 'CORE_fittedparameters_percmodel1240_respmodel15_fractrain0_20190101T162309.mat'; %'GBM_config_40BO_al4_alvar4'
% 'CORE_fittedparameters_percmodel1245_respmodel15_fractrain0_20190101T162309.mat'; %'GBM_config_45BO_al4_alvar4'
% 'CORE_fittedparameters_percmodel1250_respmodel15_fractrain0_20190101T162309.mat'; %'GBM_config_50BO_al4_alvar4' 

% empirical priors
%'CORE_fittedparameters_percmodel12_respmodel15_fractrain1_20190106T120536_it5.mat'; %'GBM_config_BO_al4_alvar1'
%'CORE_fittedparameters_percmodel12_respmodel15_fractrain1_20190108T132736_it5.mat'; %'GBM_config_BO_al4_alvar2'
%'CORE_fittedparameters_percmodel12_respmodel15_fractrain1_20190104T200202_it7.mat'; %'GBM_config_BO_al4_alvar4'
%'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190111T203212_it2.mat'; % 'GBM_config_CORE_percmodel3_BOpriors'
%'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190111T203212_it3.mat'; % 'GBM_config_CORE_percmodel3_BOpriors'
%'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190111T203212_it4.mat'; % 'GBM_config_CORE_percmodel3_BOpriors'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it1.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it2.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it3.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it4.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it5.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it6.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it7.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it8.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it9.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it10.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it11.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it12.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it13.mat'
'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it14.mat'
};

ranksum_all_p = struct;
ranksum_all_z = struct;
medians_all = struct;
mva_all = struct;
for f=1:length(fnames)

    load(fullfile(S.path.hgf,'fitted',fnames{f}));
    
    % which trial-wise stats to output for traj? Averages over trials
    % without weighting by condition
    S.summary_stats = {};%{'trial_mean','trial_std','trial_absmean','trial_absstd'}; % for all traj

    % which traj to outputs to mean over conditions? First averages for
    % each condition, before creating further averages between conditions
    %S.condmean = {'PL_muhat_1','PL_dau','PL_sahat_1','PL_muhat_2','PL_sahat_2','PL_muhat_3','PL_sahat_3'};
    %S.condmean = {'PL_dau','PL_da_1','PL_da_2','PL_epsi_1','PL_epsi_2','PL_epsi_3','PL_psi_1','PL_psi_2','PL_psi_3','PL_sa_1','PL_sa_2','PL_sa_3','PL_sahat_1','PL_sahat_2','PL_sahat_3','PL_mu_1','PL_mu_2','PL_mu_3','PL_muhat_1','PL_muhat_2','PL_muhat_3'};
    S.condmean = {'PL_dau','PL_da_1','PL_da_2','PL_epsi_1','PL_epsi_2','PL_epsi_3'};
    %S.condmean = {'PL_dau','PL_da_1','PL_da_2'};
    
%% for trajectory plots (means and variances on separate rows of S.condmean) 
%     S.condmean = {'PL_dau','PL_da_1','PL_da_2'; % for plotting traj: means
%             'PL_psi_1','PL_psi_2','PL_psi_3'}; % for plotting traj: precisions
%     S.condmean = {'PL_mu_1','PL_mu_2','PL_mu_3';
%         'PL_sa_1','PL_sa_2','PL_sa_3'};

%%
    %S.condmean = {};

    % event numbers for each condition (1st = change, 2nd = no change)
%     S.cond.CP10_Sideaff_DC1 = [1 3]; %0.1 prob, 1 digit change, left hand
%     S.cond.CP10_Sideaff_DC3 = [2 4]; %0.1 prob, 3 digit change, left hand
%     S.cond.CP10_Sideunaff_DC1 = [5 7]; %0.1 prob, 1 digit change, right hand
%     S.cond.CP10_Sideunaff_DC3 = [6 8]; %0.1 prob, 3 digit change, right hand
%     S.cond.CP30_Sideaff_DC1 = [9 11]; %0.3 prob, 1 digit change, left hand
%     S.cond.CP30_Sideaff_DC3 = [10 12]; %0.3 prob, 3 digit change, left hand
%     S.cond.CP30_Sideunaff_DC1 = [13 15]; %0.3 prob, 1 digit change, right hand
%     S.cond.CP30_Sideunaff_DC3 = [14 16]; %0.3 prob, 3 digit change, right hand
%     S.cond.CP50_Sideaff_DC1 = [17 19]; %0.5 prob, 1 digit change, left hand
%     S.cond.CP50_Sideaff_DC3 = [18 20]; %0.5 prob, 3 digit change, left hand
%     S.cond.CP50_Sideunaff_DC1 = [21 23]; %0.5 prob, 1 digit change, right hand
%     S.cond.CP50_Sideunaff_DC3 = [22 24]; %0.5 prob, 3 digit change, right hand

%% marginal means of interest
%      S.cond.condmean = [1:24]; %0.5 prob, 3 digit change, right hand
%     S.cond.odd = [1 5 9 13 17 21 2 6 10 14 18 22];  
%     S.cond.stan = [3 7 11 15 19 23 4 8 12 16 20 24];  
%      S.cond.odd_stan = {[1 5 9 13 17 21 2 6 10 14 18 22],[3 7 11 15 19 23 4 8 12 16 20 24]}; % expecting grp effects to depend on oddball vs. standard
%      S.cond.odd_CP10 = [1:2 5:6]; 
%      S.cond.odd_CP50 = [17:18 21:22]; 
%      S.cond.stan_CP10 = [3:4 7:8]; 
%      S.cond.stan_CP50 = [19:20, 23:24]; 
%      S.cond.odd_CP = {[1:2 5:6],[17:18 21:22]}; % larger oddball effects for CP10 than CP50. Grp effects expected.
%      S.cond.odd_DC = {[1 5 9 13 17 21],[2 6 10 14 18 22]}; % small DC effects on odd vs. stan 
%      S.cond.stan_CP = {[3:4 7:8],[19:20, 23:24]}; % larger oddball effects for CP10 than CP50. Grp effects expected.
%      S.cond.stan_DC = {[3 7 11 15 19 23],[4 8 12 16 20 24]}; % small DC effects on odd vs. stan 
%      S.cond.odd_stan_CP = {[1:2 5:6],[3:4 7:8],[17:18 21:22],[19:20, 23:24]}; % larger oddball effects for CP10 than CP50. Grp effects expected.
%      S.cond.odd_stan_CP10 = {[1:2 5:6],[3:4 7:8]}; % larger oddball effects for CP10 than CP50. Grp effects expected.
%      S.cond.odd_stan_CP50 = {[17:18 21:22],[19:20, 23:24]}; % larger oddball effects for CP10 than CP50. Grp effects expected.
%      S.cond.odd_stan_DC = {[1 5 9 13 17 21],[3 7 11 15 19 23],[2 6 10 14 18 22],[4 8 12 16 20 24]}; % small DC effects on odd vs. stan 

%% marginal means for violin plots
    S.cond.CP10_Odd_DC1 = [1 5];
    S.cond.CP10_Odd_DC3 = [2 6];
    S.cond.CP30_Odd_DC1 = [9 13];
    S.cond.CP30_Odd_DC3 = [10 14];
    S.cond.CP50_Odd_DC1 = [17 21];
    S.cond.CP50_Odd_DC3 = [18 22];
    S.colour_code = [1 2 1 2 1 2];
    S.x_pos = [1.2 1.8 3.2 3.8 5.2 5.8]; % x position of each violin
    S.abs = 0;

%% each cell of design: for outputting to excel for SPSS, or for plotting
%       S.cond.CP10_SideA_Odd_DC1 = 1;
%       S.cond.CP10_SideA_Odd_DC3 = 2;
%       S.cond.CP10_SideA_Stan_DC1 = 3;
%       S.cond.CP10_SideA_Stan_DC3 = 4;
%       S.cond.CP10_SideU_Odd_DC1 = 5;
%       S.cond.CP10_SideU_Odd_DC3 = 6;
%       S.cond.CP10_SideU_Stan_DC1 = 7;
%       S.cond.CP10_SideU_Stan_DC3 = 8;
%       S.cond.CP30_SideA_Odd_DC1 = 9;
%       S.cond.CP30_SideA_Odd_DC3 = 10;
%       S.cond.CP30_SideA_Stan_DC1 = 11;
%       S.cond.CP30_SideA_Stan_DC3 = 12;
%       S.cond.CP30_SideU_Odd_DC1 = 13;
%       S.cond.CP30_SideU_Odd_DC3 = 14;
%       S.cond.CP30_SideU_Stan_DC1 = 15;
%       S.cond.CP30_SideU_Stan_DC3 = 16;
%       S.cond.CP50_SideA_Odd_DC1 = 17;
%       S.cond.CP50_SideA_Odd_DC3 = 18;
%       S.cond.CP50_SideA_Stan_DC1 = 19;
%       S.cond.CP50_SideA_Stan_DC3 = 20;
%       S.cond.CP50_SideU_Odd_DC1 = 21;
%       S.cond.CP50_SideU_Odd_DC3 = 22;
%       S.cond.CP50_SideU_Stan_DC1 = 23;
%       S.cond.CP50_SideU_Stan_DC3 = 24;
%     S.colour_code = repmat([1 2],1,12);
%     S.x_pos = [1:8 13:20 25:32]; % x position of each violin
%     S.abs = 0;

%%     
%     S.cond.condmean_odd_stan = {[1 2 5 6 9 10 13 14 17 18 21 22],[3 4 7 8 11 12 15 16 19 20 23 24]}; %0.5 prob, 3 digit change, right hand
%     S.cond.DC1_odd_stan = {[1 5 9 13 17 21],[3 7 11 15 19 23]}; %0.5 prob, 3 digit change, right hand
%     S.cond.DC3_odd_stan = {[2 6 10 14 18 22],[4 8 12 16 20 24]}; %0.5 prob, 3 digit change, right hand
%     S.cond.CP10_odd_stan = {[1 2 5 6],[3 4 7 8]}; %0.5 prob, 3 digit change, right hand
%     S.cond.CP30_odd_stan = {[9 10 13 14],[11 12 15 16]}; %0.5 prob, 3 digit change, right hand
%     S.cond.CP50_odd_stan = {[17 18 21 22],[19 20 23 24]}; %0.5 prob, 3 digit change, right hand
%     S.cond.sideaff_odd_stan = {[1 2 9 10 17 18],[3 4 11 12 19 20]}; %0.5 prob, 3 digit change, right hand
%     S.cond.sideunaff_odd_stan = {[5 6 13 14 21 22],[7 8 15 16 23 24]}; %0.5 prob, 3 digit change, right hand
    
    % 
    % S.cond.CP10_Sideaff_DC1_ch = [1]; %0.1 prob, 1 digit change, left hand
    % S.cond.CP10_Sideaff_DC3_ch = [2]; %0.1 prob, 3 digit change, left hand
    % S.cond.CP10_Sideunaff_DC1_ch = [5]; %0.1 prob, 1 digit change, right hand
    % S.cond.CP10_Sideunaff_DC3_ch = [6]; %0.1 prob, 3 digit change, right hand
    % S.cond.CP30_Sideaff_DC1_ch = [9]; %0.3 prob, 1 digit change, left hand
    % S.cond.CP30_Sideaff_DC3_ch = [10]; %0.3 prob, 3 digit change, left hand
    % S.cond.CP30_Sideunaff_DC1_ch = [13]; %0.3 prob, 1 digit change, right hand
    % S.cond.CP30_Sideunaff_DC3_ch = [14]; %0.3 prob, 3 digit change, right hand
    % S.cond.CP50_Sideaff_DC1_ch = [17]; %0.5 prob, 1 digit change, left hand
    % S.cond.CP50_Sideaff_DC3_ch = [18]; %0.5 prob, 3 digit change, left hand
    % S.cond.CP50_Sideunaff_DC1_ch = [21]; %0.5 prob, 1 digit change, right hand
    % S.cond.CP50_Sideunaff_DC3_ch = [22]; %0.5 prob, 3 digit change, right hand
    % 
    % S.cond.CP10_Sideaff_DC1_noch = [3]; %0.1 prob, 1 digit change, left hand
    % S.cond.CP10_Sideaff_DC3_noch = [4]; %0.1 prob, 3 digit change, left hand
    % S.cond.CP10_Sideunaff_DC1_noch = [7]; %0.1 prob, 1 digit change, right hand
    % S.cond.CP10_Sideunaff_DC3_noch = [8]; %0.1 prob, 3 digit change, right hand
    % S.cond.CP30_Sideaff_DC1_noch = [11]; %0.3 prob, 1 digit change, left hand
    % S.cond.CP30_Sideaff_DC3_noch = [12]; %0.3 prob, 3 digit change, left hand
    % S.cond.CP30_Sideunaff_DC1_noch = [15]; %0.3 prob, 1 digit change, right hand
    % S.cond.CP30_Sideunaff_DC3_noch = [16]; %0.3 prob, 3 digit change, right hand
    % S.cond.CP50_Sideaff_DC1_noch = [19]; %0.5 prob, 1 digit change, left hand
    % S.cond.CP50_Sideaff_DC3_noch = [20]; %0.5 prob, 3 digit change, left hand
    % S.cond.CP50_Sideunaff_DC1_noch = [23]; %0.5 prob, 1 digit change, right hand
    % S.cond.CP50_Sideunaff_DC3_noch = [24]; %0.5 prob, 3 digit change, right hand


    %% run stats
    if ~exist('D_fit','var') && exist('D_sim','var')
        D_fit=D_sim;
        D_fit.HGF.fit = D_fit.HGF.sim;
    end
    [out.T,out.traj,out.param,out.rt] = CORE_extract_HGF_results(D_fit,S);
    
    %% save table
    if 0 
        sname = strrep(fullfile(S.path.hgf,'fitted',fnames{f}),'.mat','_table.xlsx');
        writetable(out.T,sname)
    end
    
    if length(D_fit)>1
        S.lda = {'traj_conds'};%{'para','traj_trials','traj_conds'};
        S.nperm=0;
        %[out.stats] = CORE_HGF_groupstatistics(out.T,{out.traj},S);
        [out.stats] = CORE_HGF_groupstatistics(out.T,{out.traj,S.condmean},S);
        %[out.stats] = CORE_HGF_groupstatistics(out.T,{},S);
    end
    
    %% combined table
    if isempty(fieldnames(ranksum_all_p))
        ranksum_all_p = out.stats.ranksum.p;
    else
        ranksum_all_p = [ranksum_all_p,out.stats.ranksum.p];
    end
    
    if isempty(fieldnames(ranksum_all_z))
        ranksum_all_z = out.stats.ranksum.z;
    else
        ranksum_all_z = [ranksum_all_z,out.stats.ranksum.z];
    end
    
    %% combined table
    try
    if isempty(fieldnames(mva_all))
        mva_all = out.stats.mvc.traj.error;
    else
        mva_all = [mva_all,out.stats.mvc.traj.error];
    end
    end
    
    %% group medians
    [grps,~,grpind] = unique(out.T.groups);
    hdr = out.T.Properties.VariableNames;
    out.median = table;
    for i = 1:size(out.T,2)
        if isnumeric(out.T.(hdr{i})(1))
            for g = 1:length(grps)
                out.median.(hdr{i})(g) = median(out.T.(hdr{i})(grpind==g));
            end
        end
    end
    %% combined table
    if isempty(fieldnames(medians_all))
        medians_all = table2struct(out.median);
    else
        medians_all = [medians_all,table2struct(out.median)];
    end

    %% save
    if 0
        save(fullfile(S.path.hgf,strrep(fnames{f},'fittedparameters','groupstatistics')), 'out');
        sname=strrep(strrep(fname,'fittedparameters','statstable'),'.mat','.xlsx');
        writetable(out.T,fullfile(S.path.hgf,sname));
    end

    %% plot group differences for selected variables
    if 0
        close all
        %vs={'PL_dau_odd','PL_dau_stan','PL_dau_odd_stan','PL_dau_odd_CP','PL_dau_stan_CP','PL_dau_odd_stan_CP','PL_epsi_1_odd','PL_epsi_2_odd','PL_epsi_3_odd','PL_epsi_1_stan','PL_epsi_2_stan','PL_epsi_3_stan','PL_psi_1_condmean','PL_psi_2_condmean','PL_psi_3_condmean','PL_psi_1_odd','PL_psi_2_odd','PL_psi_3_odd','PL_psi_1_stan','PL_psi_2_stan','PL_psi_3_stan','PL_epsi_1_odd_CP','PL_epsi_1_stan_CP','PL_epsi_1_odd_stan_CP','PL_dau_odd_DC','PL_dau_stan_DC','PL_dau_odd_stan_DC','PL_da_1_odd','PL_da_1_odd_CP','PL_da_1_odd_DC','PL_epsi_2_odd_CP','PL_epsi_2_stan_CP','PL_epsi_2_odd_stan_CP','PL_da_2_odd','PL_da_2_odd_CP','PL_da_2_odd_DC','PL_epsi_3_odd_CP','PL_epsi_3_stan_CP','PL_epsi_3_odd_stan_CP'};
        %vs={'like_al0_1','like_al0_2','like_al0_3','like_al0_4','PL_om_2','PL_om_3','PL_dau_mean','PL_da_1_mean','PL_da_2_mean','PL_epsi_1_mean','PL_epsi_2_mean','PL_epsi_3_mean','PL_psi_1_mean','PL_psi_2_mean','PL_psi_3_mean','PL_dau_odd','PL_dau_stan','PL_dau_odd_stan','PL_dau_odd_CP','PL_dau_stan_CP','PL_dau_odd_stan_CP','PL_dau_odd_DC','PL_dau_stan_DC','PL_dau_odd_stan_DC','PL_da_1_odd','PL_da_1_odd_CP','PL_da_1_odd_DC','PL_da_2_odd','PL_da_2_odd_CP','PL_da_2_odd_DC'};
        %vs={'like_al0','PL_om_2','PL_om_3','PL_sa_2_mean','PL_ud_1_mean','PL_ud_2_mean','PL_dau_mean','PL_epsi_2_mean','PL_epsi_3_mean'};
        vs={'PL_dau_condmean','PL_da_1_condmean','PL_da_2_condmean','PL_psi_1_condmean','PL_psi_2_condmean','PL_psi_3_condmean'};
        for v=1:length(vs)
            figure
            %scatter(grpind,out.T.(vs{v}))
            for g = 1:length(grps)
                X{g}=out.T.(vs{v})(grpind==g);
            end
            n_rainclouds(X,cb(1:g,:))
            %set(gca,'XTick',1:2)
            %set(gca,'XTickLabel',grps)
            title([vs{v} ', p = ' num2str(out.stats.ranksum.p.(vs{v}))])
        end
    end
    
    %% HGF trajectories
    if 0
        close all
        subplot_on = 1;
        plot_variance_type = 'mean_of_traj'; %{'mean_of_traj','var_over_subs'};
        var_pi = 'pi'; % precision (pi) or variance (var)?
        
        % group indices
        for g = 1:length(grps)
            % Set up display
            scrsz = get(0,'screenSize');
            outerpos = [0.2*scrsz(3),0.2*scrsz(4),0.8*scrsz(3),0.8*scrsz(4)];
            figure(...
                'OuterPosition', outerpos,...
                'Name', ['HGF trajectories, ' grps{g}]);

            vs=S.condmean(1,:);
            [cb] = cbrewer('qual', 'Set1', length(vs), 'pchip');

            leg={''};
            for v=1:length(vs)

                if subplot_on
                    ax(g,v)=subplot(length(vs),1,length(vs)-v+1);
                else
                    hold on
                end

                % get data
                dat = cat(2,out.traj(:).(vs{v}));
                mn = mean(dat(:,grpind==g),2)';
                ns = length(grpind==g);
                ts=1:length(mn);

                if strcmp(plot_variance_type,'var_over_subs')
                    % get spread
                    sd = std(dat(:,grpind==g),[],2)';
                    sem = sd/sqrt(ns);              % Standard Error
                    tscore = tinv(0.05,ns-1);      % T-Score
                    CI = tscore*sem;                % Confidence Intervals
                    % choose spread type
                    spr = sd/2;
                elseif strcmp(plot_variance_type,'mean_of_traj')
                    vs_var=S.condmean{2,v};
                    % get data
                    dat_var = cat(2,out.traj(:).(vs_var));
                    dat_var = dat_var(:,grpind==g);
                    if strcmp(var_pi,'pi')
                        dat_var = 1./dat_var;
                    end
                    spr = mean(dat_var,2)'/2;
                end

                upper = mn+spr;
                lower = mn-spr;

                % plot spread
                fill([ts, fliplr(ts)], [(upper), fliplr((lower))], ...
                cb(v,:), 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
                hold on

                % plot mean
                ln(v) = plot(ts,mn,'Color',cb(v,:));
                lim(g,v,:)=axis;
                if subplot_on
                    ylabel(vs{v})
                    if v==1
                        xlabel('Trial number')
                    end
                end
                hold off
            end
            % plot inputs/design
            %plot(ts, D_fit(1).HGF.u(:,1), '.', 'Color', [0 0.6 0]); % inputs
            if ~subplot_on
                legend(ln,vs)
                ylabel('Prediction error')
                xlabel('Trial number')
            end
        end
        % equate axes between figures
        maxlim(:,[2,4]) = squeeze(max(lim(:,:,[2 4]),[],1));
        maxlim(:,[1,3]) = squeeze(min(lim(:,:,[1 3]),[],1));
        for g = 1:length(grps)
            for v=1:length(vs)
                axis(ax(g,v),maxlim(v,:));
            end
        end
    end
    


    %% condition means of HGF trajectories
    if 1
        vs=S.condmean;
        ucol = unique(S.colour_code);
        [cb] = cbrewer('qual', 'Set1', length(ucol), 'pchip');
        for v=1:length(vs)

            % get column indices
            colnames = out.T.Properties.VariableNames;
            vss{v} = colnames(~cellfun(@isempty,strfind(colnames,vs{v})));
            
            % get data
            dat = out.T(:,vss{v});
            xlab = dat.Properties.VariableNames;
            if S.abs
                plotdat = abs(dat{:,:});
            else
                plotdat = dat{:,:};
            end

            % plot
            figure('Name','All groups')
            %bar(xlab,dat{1,:})
            violinplot(plotdat,xlab,S.x_pos,'ViolinColor',cb(S.colour_code,:),'ViolinAlpha',0.6,'EdgeColor',[0.8 0.8 0.8],'MedianColor',[0 0 0])
            title(vs{v})
            
            % separate plots by group indices
            for g = 1:length(grps)

                % get data
                dat = out.T(grpind==g,vss{v});
                xlab = dat.Properties.VariableNames;
                if S.abs
                    plotdat = abs(dat{:,:});
                else
                    plotdat = dat{:,:};
                end

                % plot
                figure('Name',grps{g})
                %bar(xlab,dat{1,:})
                violinplot(plotdat,xlab,S.x_pos,'ViolinColor',cb(S.colour_code,:),'ViolinAlpha',0.6,'EdgeColor',[0.8 0.8 0.8],'MedianColor',[0 0 0])
                title(vs{v})
            end
        end
    end
end
try [mva_all(:).fname] = deal(fnames{:});end
[ranksum_all_p(:).fname] = deal(fnames{:});

% flag significance
fd=fieldnames(ranksum_all_p);
ii = length(ranksum_all_p)+1;
for fn = 1:length(fd)
    ranksum_all_p(ii).(fd{fn}) = sum([ranksum_all_p(1:ii-1).(fd{fn})]<0.05)>length(ranksum_all_p)/10;
end
