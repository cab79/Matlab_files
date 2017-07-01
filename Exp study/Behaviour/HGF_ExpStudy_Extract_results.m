cd('M:\Matlab\ExpStudy');
clear all

max_stim = 360;
rname = 'ExpStudy_HGF_results_v2.mat';
load(rname);

subs = fieldnames(results);
nsub =length(subs);

al2 = nan(nsub,1);
alpf2 = nan(nsub,1);
om1_1 = nan(nsub,1);
om1_2 = nan(nsub,1);
om1_3 = nan(nsub,1);
om2_1 = nan(nsub,1);
om2_2 = nan(nsub,1);
om2_3 = nan(nsub,1);
dau = nan(nsub,2);
da1 = nan(nsub,2);
da2 = nan(nsub,2);
da3 = nan(nsub,2);
ud1 = nan(nsub,2,6);
ud2 = nan(nsub,2,6);
ud3 = nan(nsub,2,6);
wt1 = nan(nsub,2);
wt2 = nan(nsub,2);
wt3 = nan(nsub,2);
psi1 = nan(nsub,2);
psi2 = nan(nsub,2);
psi3 = nan(nsub,2);
epsi1 = nan(nsub,2);
epsi2 = nan(nsub,2);
epsi3 = nan(nsub,2);
ud_traj = nan(nsub,2,max_stim);

for s = 1:nsub;
    sub = subs{s};
    stimlen = length(results.(sub).traj.dau);
    cond = results.(sub).cond;
    al2(s,1) = results.(sub).p_prc.al2;
    alpf2(s,1) = results.(sub).p_prc.alpf2;
    om1_1(s,1) = results.(sub).p_prc.om1(1);
    om1_2(s,1) = results.(sub).p_prc.om1(2);
    om1_3(s,1) = results.(sub).p_prc.om1(3);
    om2_1(s,1) = results.(sub).p_prc.om2(1);
    om2_2(s,1) = results.(sub).p_prc.om2(2);
    om2_3(s,1) = results.(sub).p_prc.om2(3);
    for i = 1:2
        ind = i:2:(stimlen-2+i);
        dau(s,i) = nanmean(results.(sub).traj.dau(ind));
        da1(s,i) = nanmean(results.(sub).traj.da(ind,1),1);
        da2(s,i) = nanmean(results.(sub).traj.da(ind,2),1);
        da3(s,i) = nanmean(results.(sub).traj.da(ind,3),1);
        for c = 1:6
            cind = find(cond == c)';
            ud1(s,i,c) = nanmean(results.(sub).traj.ud(ind(cind),1),1);
            ud2(s,i,c) = nanmean(results.(sub).traj.ud(ind(cind),2),1);
            ud3(s,i,c) = nanmean(results.(sub).traj.ud(ind(cind),3),1);
        end
        ud_traj(s,i,1:stimlen/2) = results.(sub).traj.ud(ind,1);
        wt1(s,i) = nanmean(results.(sub).traj.wt(ind,1),1);
        wt2(s,i) = nanmean(results.(sub).traj.wt(ind,2),1);
        wt3(s,i) = nanmean(results.(sub).traj.wt(ind,3),1);
        psi1(s,i) = nanmean(results.(sub).traj.psi(ind,1),1);
        psi2(s,i) = nanmean(results.(sub).traj.psi(ind,2),1);
        psi3(s,i) = nanmean(results.(sub).traj.psi(ind,3),1);
        epsi1(s,i) = nanmean(results.(sub).traj.epsi(ind,1),1);
        epsi2(s,i) = nanmean(results.(sub).traj.epsi(ind,2),1);
        epsi3(s,i) = nanmean(results.(sub).traj.epsi(ind,3),1);
    end
end

dau1 = dau(:,1);
dau2 = dau(:,2);
da1_1 = da1(:,1);
da2_1 = da1(:,2);
da1_2 = da2(:,1);
da2_2 = da2(:,2);
da1_3 = da3(:,1);
da2_3 = da3(:,2);

ud1_1_1 = ud1(:,1,1);
ud2_1_1 = ud1(:,2,1);
ud1_2_1 = ud2(:,1,1);
ud2_2_1 = ud2(:,2,1);
ud1_3_1 = ud3(:,1,1);
ud2_3_1 = ud3(:,2,1);
ud1_1_2 = ud1(:,1,2);
ud2_1_2 = ud1(:,2,2);
ud1_2_2 = ud2(:,1,2);
ud2_2_2 = ud2(:,2,2);
ud1_3_2 = ud3(:,1,2);
ud2_3_2 = ud3(:,2,2);
ud1_1_3 = ud1(:,1,3);
ud2_1_3 = ud1(:,2,3);
ud1_2_3 = ud2(:,1,3);
ud2_2_3 = ud2(:,2,3);
ud1_3_3 = ud3(:,1,3);
ud2_3_3 = ud3(:,2,3);
ud1_1_4 = ud1(:,1,4);
ud2_1_4 = ud1(:,2,4);
ud1_2_4 = ud2(:,1,4);
ud2_2_4 = ud2(:,2,4);
ud1_3_4 = ud3(:,1,4);
ud2_3_4 = ud3(:,2,4);
ud1_1_5 = ud1(:,1,5);
ud2_1_5 = ud1(:,2,5);
ud1_2_5 = ud2(:,1,5);
ud2_2_5 = ud2(:,2,5);
ud1_3_5 = ud3(:,1,5);
ud2_3_5 = ud3(:,2,5);
ud1_1_6 = ud1(:,1,6);
ud2_1_6 = ud1(:,2,6);
ud1_2_6 = ud2(:,1,6);
ud2_2_6 = ud2(:,2,6);
ud1_3_6 = ud3(:,1,6);
ud2_3_6 = ud3(:,2,6);

wt1_1 = wt1(:,1);
wt2_1 = wt1(:,2);
wt1_2 = wt2(:,1);
wt2_2 = wt2(:,2);
wt1_3 = wt3(:,1);
wt2_3 = wt3(:,2);
psi1_1 = psi1(:,1);
psi2_1 = psi1(:,2);
psi1_2 = psi2(:,1);
psi2_2 = psi2(:,2);
psi1_3 = psi3(:,1);
psi2_3 = psi3(:,2);
epsi1_1 = epsi1(:,1);
epsi2_1 = epsi1(:,2);
epsi1_2 = epsi2(:,1);
epsi2_2 = epsi2(:,2);
epsi1_3 = epsi3(:,1);
epsi2_3 = epsi3(:,2);
        
T1 = table(al2,alpf2,om1_1,om1_2,om1_3,om2_1,om2_2,om2_3,dau1,dau2,da1_1,da1_2,da1_3,da2_1,da2_2,da2_3,wt1_1,wt1_2,wt1_3,wt2_1,wt2_2,wt2_3,psi1_1,psi1_2,psi1_3,psi2_1,psi2_2,psi2_3,epsi1_1,epsi1_2,epsi1_3,epsi2_1,epsi2_2,epsi2_3,'RowNames',subs);
T2 = table(ud1_1_1,ud1_1_2,ud1_1_3,ud1_1_4,ud1_1_5,ud1_1_6,ud1_2_1,ud1_2_2,ud1_2_3,ud1_2_4,ud1_2_5,ud1_2_6,ud1_3_1,ud1_3_2,ud1_3_3,ud1_3_4,ud1_3_5,ud1_3_6,ud2_1_1,ud2_1_2,ud2_1_3,ud2_1_4,ud2_1_5,ud2_1_6,ud2_2_1,ud2_2_2,ud2_2_3,ud2_2_4,ud2_2_5,ud2_2_6,ud2_3_1,ud2_3_2,ud2_3_3,ud2_3_4,ud2_3_5,ud2_3_6,'RowNames',subs);

tname = 'ExpStudy_HGF_results_v3.xlsx';
writetable(T1,tname);

tname = 'ExpStudy_HGF_results_v3_ud.xlsx';
writetable(T2,tname);


