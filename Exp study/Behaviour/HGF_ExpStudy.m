clear all
pth=('M:\Matlab\ExpStudy');
cd(pth);
ci = 3;
pro_ana = 1; % proportion of data to analyse

results = struct;

ncolhead=3;
[data hdr raw]=xlsread('Pain_ratings_HGF.xls');
rname = 'ExpStudy_HGF_results_v4.mat';
if exist(rname,'file');
    load(rname);
end

nsub = length(hdr)-ncolhead;

intall = data(:,find(ismember(hdr,'int'))); 
uncertall = data(:,find(ismember(hdr,'uncert')));
condall = data(:,find(ismember(hdr,'cond')));

for s=floor((nsub/3)*(ci-1)+1):nsub
    %sort([1:16],'descend')
    
    
    % load inputs
    subname = hdr{ncolhead+s}
    prn=data(:,ncolhead+s);
    
    if isfield(results,subname)
        continue
    end
    
    % find first and last non-NaN
    pr_idx0 = find(~isnan(prn));
    if ~isempty(pr_idx0); 
        pr_idx = pr_idx0([1,end]);
    else
        pr_idx = [1,length(prn)];
    end
    
    % data for analysis
    pr = prn(pr_idx(1):pr_idx(2));
    int = intall(pr_idx(1):pr_idx(2));
    uncert = uncertall(pr_idx(1):pr_idx(2));
    cond = condall(pr_idx(1):pr_idx(2));
    
    % remove trials from extended periods of no data
    prin = isnan(pr);
    prin_idx=zeros(length(prin),1);
    for i = 3:length(prin)-3
        if (prin(i) && prin(i+1) && prin(i+2)) || prin(i) && prin(i-1) && prin(i-2)
            prin_idx(i,1)=1;
        end
    end
    pr(find(prin_idx))=[];
    int(find(prin_idx))=[];
    uncert(find(prin_idx))=[];
    cond(find(prin_idx))=[];
        
    % For each certain condition, use RW model to generate trajectory of cue values over time 
    cueval = nan(size(pr));
    for c = 4:6
        cond_idx = find(cond==c);
        pr_cond = pr(cond_idx);
        u=pr_cond;
        [traj] = rw_model(u);
        cueval(cond_idx)=traj.v;
    end
    
    % For each uncertain condition, use ema to generate trajectory of nociception over time 
    nocval = nan(size(pr));
    for in = [3 5 7]
        int_idx = find(int==in);
        pr_int = pr(int_idx);
        pr_av = nan(size(pr_int));
        for i = 1:length(pr_int)
             if isnan(pr_int(i)) && i==1
                 pr_int(i) = nanmean(pr_int(i:i+5));
             elseif isnan(pr_int(i)) && i>1
                 pr_int(i) = nanmean(pr_int(max(1,i-5):i));
             end
             ema_all = ema(pr_int(1:i),i);
             pr_av(i) = ema_all(end);
        end
        nocval(int_idx)=pr_av;
    end

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
    
    data_ana = floor(length(u)*pro_ana);

    % prc: perceptual; obs:observation; opt:optimisation
    prc_model = 'tapas_hgf_salience_2input_config_2level_v4';
    obs_model = 'tapas_gaussian_obs_post_config';
    %obs_model = 'tapas_bayes_optimal_config'; % ignores responses
    opt_algo = 'tapas_quasinewton_optim_config';
    bopars = tapas_fitModel(y(1:data_ana,:), u(1:data_ana,:), prc_model, obs_model, opt_algo);
    %tapas_hgf_salience_2input_plotTraj(bopars)
    
    if exist(rname,'file'); load(rname); end;
    
    results.(subname)=bopars;
    results.(subname).u=u;
    results.(subname).y=y;
    results.(subname).pr=pr;
    results.(subname).cond=cond;
    results.(subname).int=int;
    results.(subname).uncert=uncert;
    
    save(rname,'results');

end


