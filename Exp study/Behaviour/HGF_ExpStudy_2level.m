clear all

%for r = 1:3

rawpth=('M:\Matlab\ExpStudy\Behaviour\Raw');
pth=('M:\Matlab\ExpStudy\Behaviour\Results');
condpth=('M:\Matlab\ExpStudy\Behaviour\Results\cond numbers');
cd(pth);
ci = 1;
smoothpara =1;

results = struct;
rname = ['ExpStudy_HGF_results_LLS' num2str(smoothpara) '.mat'];

files = dir(fullfile(condpth,'*conds.mat'));
nsub = length(files);

for s=sort(floor((nsub/3)*(ci-1)+1):nsub, 'ascend')
    %sort([1:16],'descend')
    
    % load inputs
    fname = files(s).name;
    C = strsplit(fname,'_');
    subname = C{1};
    load(fullfile(condpth,fname));
    if exist(rname,'file'); load(rname); end;

    
    %nocval1 = nan(size(pr));
    %for c=[1:3]
    %    cond_idx = find(cond==c);
    %    cond2_idx = find(cond==c+3);
    %    int_idx = find(int==pint(c));
    %    nocval1(cond_idx) = pr(cond_idx);
    %    pr_int = nocval1(int_idx);
    %    for i = 1:length(pr_int)
    %        if i==1 && isnan(pr_int(i))
    %            pr_int(i)=pint(cond(int_idx(i)));
    %        elseif isnan(pr_int(i))
   %            pr_int(i)=pr_int(i-1);
   %         end
   %     end
   %     pr_sm1 = smooth(pr_int,5,'lowess');
   %     %pr_sm2 = smooth(pr_int,5,'moving');
   %     %pr_sm3 = smooth(pr_int,80,'lowess');
   %     %pr_sm4 = smooth(pr_int,10,'moving');
   %     nocval1(int_idx)=pr_sm1;
   % end
   
    exp_effect = nan(3,1);
    for c = 1:3
        condC_idx = find(cond==c+3);
        condU_idx = find(cond==c);
        exp_effect(c,1) = (nanmean(pr(condC_idx))-nanmean(pr(condU_idx)))*100/nanmean(pr(condU_idx));
    end
   
   % For each certain condition, use RW model to generate trajectory of cue values over time 
 %   cueval = nan(size(pr));
 %   for c = 4:6
 %       cond_idx = find(cond==c);
 %       pr_cond = pr(cond_idx);
 %       [traj] = rw_model(pr_cond);
 %       cueval(cond_idx(1))=pr_cond(1);
 %       cueval(cond_idx(2:end))=traj.v(1:end-1); % cueval is the value of the previous cue
      
        %pr_cond = nocval(cond_idx);
        %[traj] = rw_model(pr_cond);
        %cueval(cond_idx)=traj.v;
        %plot(cueval(cond_idx),'r'); hold on
       
 %   end
    
    % calculate Cueval for uncertain conditions
  %  cond_idx = find(uncert==1);
  %  pr_cond = pr(cond_idx);
  %  [traj] = rw_model(pr_cond);
  %  cueval(cond_idx(1))=pr_cond(1);
  %  cueval(cond_idx(2:end))=traj.v(1:end-1); % cueval is the value of the previous cue

    % For each intensity condition, use ema to generate trajectory of nociception over time 
    %nocval1 = nan(size(pr));
    nocval = nan(size(pr));
    %nocval3 = nan(size(pr));
    %nocval4 = nan(size(pr));
    %nocval7 = nan(size(pr));
    for in = [3 5 7]
        int_idx = find(int==in);
        pr_int = pr(int_idx);
        %[traj] = rw_model(pr_int);
        %nocval1(int_idx)=traj.v;
        
        %pr_av = nan(size(pr_int));
        %for i = 1:length(pr_int)
        %     ema_all = ema(pr_int(1:i),i);
        %     pr_av(i) = ema_all(end);
        %end
        
        %pr_sm1 = smooth(pr_int,1,'lowess');
        %pr_sm2 = smooth(pr_int,3,'lowess');
        %pr_sm = smooth(pr_int,smoothpara,'lowess');
        %pr_sm7 = smooth(pr_int,40,'lowess');
        p = polyfit(1:length(pr_int),pr_int',1);
        pr_sm=1:length(pr_int); pr_sm=pr_sm*p(1)+p(2);
        pr_sm=pr_sm';
        
        %pr_sm4 = pr_sm7 + (pr_int - cueval(int_idx));
        
        %nocval7(int_idx)=pr_sm7;
        %nocval3(int_idx)=pr_sm3;
        %nocval4(int_idx)=pr_sm4;
        
        %eval(['nocval(int_idx)=pr_sm' num2str(r) ';']);
        nocval(int_idx)=pr_sm;
    end
    
    cueval = nan(size(pr));
    for c = 4:6
        cond_idx = find(cond==c);
        cueval(cond_idx)=nocval(cond_idx);
    end
        
    
    
    
    %plots
    %close all
    %for c = 1:6
    %    cond_idx = find(cond==c);
    %    pr_cond = pr(cond_idx);
    %    figure
    %   plot(pr_cond,'k'); hold on
    %   plot(nocval3(cond_idx),'r');
    %   plot(nocval7(cond_idx),'y');
    %   plot(cueval(cond_idx),'b');
    %end
    %for c = 1:6
    %    cond_idx = find(cond==c);
    %    pr_cond = pr(cond_idx);
    %    figure
    %    plot(pr_cond,'k'); hold on
    %    plot(nocval4(cond_idx),'r');
    %    plot(nocval7(cond_idx),'y');
    %    plot(cueval(cond_idx),'b');
    %end
    %for c = 1:3
    %    cond_idx = find(cond==c);
    %    pr_int = pr(cond_idx);
    %    figure
    %    plot(pr_int,'k'); hold on
    %    plot(nocval(cond_idx),'r');
    %end
    %for c = 4:6
    %    cond_idx = find(cond==c);
    %    figure
    %    plot(pr(cond_idx),'k'); hold on
    %    plot(nocval(cond_idx),'r'); hold on
    %    plot(cueval(cond_idx),'b');
    %end

    % input intersperses cues with stimulus
    u = nan(length(pr)*2,2);
    y = nan(length(pr)*2,1);
    u(1:2:length(u)-1,1) = cueval;
    u(2:2:length(u),1) = nocval;
    u(1:2:length(u)-1,2) = ones(length(cueval),1);
    u(2:2:length(u),2) = 2*ones(length(nocval),1);
    y(2:2:length(u)) = pr;
    
    %convert u and y to values between 0 and 1
    u(:,1) = u(:,1)/10;
    y = y/10;

    % prc: perceptual; obs:observation; opt:optimisation
    prc_model = 'tapas_hgf_salience_2input_config_2level_v4';
    obs_model = 'tapas_gaussian_obs_post_config';
    %obs_model = 'tapas_softmax_config'; 
    opt_algo = 'tapas_quasinewton_optim_config';
    bopars = tapas_fitModel(y, u, prc_model, obs_model, opt_algo);
    %tapas_hgf_salience_2input_plotTraj(bopars)
    
    if exist(rname,'file'); load(rname); end;
    
    results.(subname)=bopars;
    results.(subname).u=u;
    results.(subname).y=y;
    results.(subname).pr=pr;
    results.(subname).cond=cond;
    results.(subname).int=int;
    results.(subname).uncert=uncert;
    results.(subname).expeffect=exp_effect;
    
    save(rname,'results');

end

%end
