clear all


Pbase = 'E:\Expectancy study\Results\Raw data';
%**************************************

cd(Pbase);
save Pbase.mat Pbase;


Dpath = spm_select(Inf, 'dir', 'Select subject directories');
Nsub = size(Dpath, 1);
save Dpath.mat Dpath;

grand_avg1=zeros(64,689);
grand_avg2=zeros(64,689);
grand_avg3=zeros(64,689);
grand_avg4=zeros(64,689);
grand_avg5=zeros(64,689);
grand_avg6=zeros(64,689);


for n = 1:Nsub
	
    cd(Pbase)
    
    load Dpath.mat;
    Dpath = deblank(Dpath(n, :));
    
    Pdata = Dpath;

    cd(Pdata);

    load('avg1_le.mat')
    load('avg2_le.mat')
    load('avg3_le.mat')
    load('avg4_le.mat')
    load('avg5_le.mat')
    load('avg6_le.mat')
    grand_avg1=grand_avg1 + mdata1;
    grand_avg2=grand_avg2 + mdata2;
    grand_avg3=grand_avg3 + mdata3;
    grand_avg4=grand_avg4 + mdata4;
    grand_avg5=grand_avg5 + mdata5;
    grand_avg6=grand_avg6 + mdata6;
    
end

cd(Pbase);
grand_avg1 = grand_avg1 ./ Nsub;
grand_avg2 = grand_avg2 ./ Nsub;
grand_avg3 = grand_avg3 ./ Nsub;
grand_avg4 = grand_avg4 ./ Nsub;
grand_avg5 = grand_avg5 ./ Nsub;
grand_avg6 = grand_avg6 ./ Nsub;
grand_U=(grand_avg3+grand_avg2+grand_avg1)/3;
grand_C=(grand_avg4+grand_avg5+grand_avg6)/3;
grand_all=(grand_avg1+grand_avg2+grand_avg3+grand_avg4+grand_avg5+grand_avg6)/6;
save grand_all_Hea.mat grand_all

t=-500:188;
plot(t, -grand_avg4(14,:), 'b', t, -grand_avg5(14,:), 'm', t, -grand_avg6(14,:), 'r', t, -grand_U(14,:), 'g');

t=-500:188;
plot(t, -grand_C(15,:), 'b',t, -grand_U(15,:), 'k');

t=-500:188;
plot(t, -grand_all(15,:), 'k')

t=-500:188;
plot(t, -mdata1(15,:), 'k')
hold on
t=-500:188;
plot(t, -grandavg(15,:), 'r')

load grand_all_FMS
t=-500:188;
plot(t, -grand_all(15,:), 'r')
hold on;
load grand_all_Hea
t=-500:188;
plot(t, -grand_all(15,:), 'b')
