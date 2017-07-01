clear all
dname=('C:\Data\CORE\Behaviour');
cd(dname);

files = dir('dt_*startblock1*');
files_ana = [1:12];

acc_all=[];
rt_all=[];
for f = files_ana

    C = strsplit(files(f).name,'_');
    
    %identify if data is split into blocks and combine
    %bfile = dir([C{1} '_' C{2} '_' C{3} '_startblock2*']);
    %if ~isempty(bfile)
        
    dt_name = files(f).name;
    RT_name = ['RT_' dt_name(4:end)];
    
    % load inputs
    %u = load('example_cont_input.txt');
    load(dt_name);
    load(RT_name);
    design=dt.design;
    u = design(3,:)'; % change = 1, no change = 0;
    comb = [u';RT(1,:)];
    
    [hand,dc,cp,bi,blockii,btypes] = blocktype(dname,dt_name);

    % event numbers for each condition (1st = change, 2nd = no change)
    cond1 = [1 3]; %0.1 prob, 1 digit change, left hand
    cond2 = [2 4]; %0.1 prob, 3 digit change, left hand
    cond3 = [5 7]; %0.1 prob, 1 digit change, right hand
    cond4 = [6 8]; %0.1 prob, 3 digit change, right hand
    cond5 = [9 11]; %0.3 prob, 1 digit change, left hand
    cond6 = [10 12]; %0.3 prob, 3 digit change, left hand
    cond7 = [13 15]; %0.3 prob, 1 digit change, right hand
    cond8 = [14 16]; %0.3 prob, 3 digit change, right hand
    cond9 = [17 19]; %0.5 prob, 1 digit change, left hand
    cond10 = [18 20]; %0.5 prob, 3 digit change, left hand
    cond11 = [21 23]; %0.5 prob, 1 digit change, right hand
    cond12 = [22 24]; %0.5 prob, 3 digit change, right hand

    ISI = 1.0;
    thresh = 0.2; % min RT that is realistic
    num_trial_lag = [1 1 1]; % number of trials over which response is allowed to lag for each of 0.1 prob, 0.3 prob, 0.5 prob.

    for i = 2:size(design,2) % ignore first trial as unlikely to incur a response
        if design(2,i)==0 % mark changes of block as being a stimulus change
            comb(1,i)=1;
        end
    end
    
    for i = 1:size(design,2)-1
        if any(cond1 == design(2,i)) || any(cond2 == design(2,i)) || any(cond3 == design(2,i)) || any(cond4 == design(2,i)) 
            lag = num_trial_lag(1);
        elseif any(cond5 == design(2,i)) || any(cond6 == design(2,i)) || any(cond7 == design(2,i)) || any(cond8 == design(2,i)) 
            lag = num_trial_lag(2);
        elseif any(cond9 == design(2,i)) || any(cond10 == design(2,i)) || any(cond11 == design(2,i)) || any(cond12 == design(2,i)) 
            lag = num_trial_lag(3);
        end
        
        comb(3,i)=0; % row 3 is the corrected response time on the corrected trial
        if comb(1,i)>0 % if the stimulus actually changed, but the response occured on the next trial, count it as occuring on the current trial but with a longer RT
            for j = 1:lag+1
                move=0;
                if comb(3,i)==0 && comb(2,i+(j-1))>0 && comb(3,i-1)~=ISI*1+comb(2,i+(j-1)) && comb(2,i+(j-1))+(j-1)*ISI>thresh
                    if comb(1,i+(j-1))==0 || j==1 % if we are considering the current trial only, or if there is no stimulus on the next trial, then it's always ok to register the response on te current trial.
                        move=1;
                    elseif comb(1,i+(j-1))>0 && j>1 %if there is a stimulus on the next trial AND
                        if comb(2,i+(j-1))>0 && comb(2,i+(j-1))<thresh % the response is less than threshold OR
                            move=1;
                        elseif (comb(1,i+j)==0 && comb(2,i+j)>0) % there is no stim on the subsequent trial but a response
                            move=1;
                        end
                    end
                end
                if move==1
                    comb(3,i)=comb(2,i+(j-1)) + (j-1)*ISI; % adds on the ISI from RT on next trial and places it on the current trial
                end
            end
        end
        if comb(1,i)==0 && comb(2,i)>0 && i>1 % if no stimulus change, but a response, 
            if ~any(comb(1,i-(j-1):i)>0) % and there was no stim change on previous trial,, so it can't be a delayed response to a stim change
                if comb(2,i)>thresh;
                    comb(3,i)=comb(2,i); % count it as a false positive response
                elseif comb(3,i-1)==0
                    comb(3,i-1)=comb(2,i)+1*ISI; % or a response on the previous trial if too fast for current trial
                end
            end
        end
    end
    
    % create labels for block transitions between hands and within hands
    u2=zeros(1,size(comb,2));
    for i = 2:length(blockii) 
        if hand(i) ~= hand(i-1) % label transitions between hands as 7
            u2(blockii(i))=7; 
        elseif dc(i) ~= dc(i-1) % label transitions within hands (from DC1 to DC3) as 5 for left, 6 for right
            if hand(i)==1
                u2(blockii(i))=5; 
            elseif hand(i)==2
                u2(blockii(i))=6;
            end
        end
    end
    comb(1,blockii)=u2(1,blockii);
    for i = 1:length(u2)
        if design(2,i)~=0
            ii = find(any(cell2mat(btypes(:,1))==design(2,i),2));
            handval = btypes{ii,2};
            dcval = btypes{ii,3};
            if handval==1 && dcval==1
                u2(i)=1;
            elseif handval==1 && dcval==2
                u2(i)=2;
            elseif handval==2 && dcval==1
                u2(i)=3;
            elseif handval==2 && dcval==2
                u2(i)=4;
            end
        end
    end
    
    % CHANGE REMAINING ZEROS IN U2 to be same as previous trial
    zi = find(u2==0);
    zi(zi==1)=[]; % ignore first trial
    u2(zi) = u2(zi-1);
    u2(1) = u2(2);

    %conds_ana = [1 3 5];
    %conds_ana = [2 4 6];
    %idx=[];
    %for i = 1:length(conds_ana)
    %    eval(['conds = cond' num2str(conds_ana(i)) ';']); 
    %    for j = 1:length(conds)
    %        idx = [idx find(design(2,:)==conds(j))];
    %    end
    %end
    %idx = sort(idx);
    %design = design(:,idx);
    %comb = comb(:,idx);  
    y = double(comb(3,:)>0); % BINARY response
    y=y';
    u=[];
    u(1,:) = comb(1,:)>0;
    %u2(u2>2) =1;
    u(2,:) = u2;
    u=u';

    %for i = 1:length(u)
    %    if u(i)==0.1
    %        u(i)=1.75;
    %    end
    %end
    
    %% analysis % correct and RTs for each condition 1-12
    transi = ismember(u2,[1:4]); %use only transitions of interest (1-4 from u2)
    % create index of cond numbers corresponding to each trial
    condi =nan(1,length(u2));
    blockii_end = [blockii length(u2)+1];
    for b = 1:length(bi)
        condi(blockii_end(b):blockii_end(b+1)-1) = bi(b);
    end
    conds = sort(unique(bi));
    pc_corr = nan(length(conds),1);
    rt_corr = nan(length(conds),1);
    for c = conds
        uyr=[];
        uyr(1,:) = u(condi==c,1);
        uyr(2,:) = y(condi==c);
        uyr(3,:) = comb(3,condi==c);
        uyr(4,:) = transi(condi==c);
        uyr = uyr(:,uyr(4,:)==1);
        pc_corr(c) = 100*sum(double(uyr(1,:)==uyr(2,:)))/size(uyr,2);
        rt_corr(c) = mean(uyr(3,intersect(find(uyr(1,:)==1), find(uyr(2,:)==1))));
    end
    acc_all(:,f)=pc_corr;
    rt_all(:,f)=rt_corr;

    %% HGF
    % prc: perceptual; obs:observation; opt:optimisation
    prc_model = 'tapas_hgf_binary_pu_CORE_config';
    %prc_model = 'tapas_hgf_binary_pu_config';
    obs_model = 'tapas_softmax_binary_CORE_config';
    %obs_model = 'tapas_bayes_optimal_binary_config'; %BAYES OPTIMAL
    opt_algo = 'tapas_quasinewton_optim_config';
   % bopars = tapas_fitModel(y, u, prc_model, obs_model, opt_algo);
    %bopars = tapas_fitModel([], u, prc_model, obs_model, opt_algo); %BAYES OPTIMAL
    %tapas_fit_plotCorr(bopars)
   % sname = [C{2} '_bopars'];
   % save(sname,'bopars');
   % tapas_hgf_binary_plotTraj(bopars)
    
    %NOTES
    % softmax doesn't like multiple columns in U
    

    %alpha = 10; %uncertain
    %alpha = 0.5; %a bit uncertain
    %alpha = 0.05; %certain
    %sim = tapas_simModel(u(2:end), 'tapas_hgf', [-0.0197 1 1 0.09 1 1 0 0 0 1 1 -1.1537 -6.2314 -6.2823 alpha], 'tapas_softmax_binary',[log(48)]);
    %tapas_hgf_plotTraj(sim)

    %signal detection values
   % input = u;
    %response=sim.y;
    %[d,beta,C] = signal_detection(input(1:end-1),response(2:end))
   % response=y;
   % [d,beta,C] = signal_detection(input,response)

    % change the responses simulated to make 0/1 changes less perceptible
    %prc_model = 'tapas_hgf_config';
    %obs_model = 'tapas_gaussian_obs_config';
    %obs_model = 'tapas_softmax_binary_config';
    %opt_algo = 'tapas_quasinewton_optim_config';
    %bopars = tapas_fitModel(sim.y, u, prc_model, obs_model, opt_algo);
    %bopars = tapas_fitModel(zeros(length(sim.y),1), u, prc_model, obs_model, opt_algo);
    %tapas_hgf_plotTraj(bopars)
end
xlswrite('accuracy.xlsx',acc_all');
xlswrite('reactiontime.xlsx',rt_all');