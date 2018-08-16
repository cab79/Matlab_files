function [S,D] = CORE_data_process(S,D)

ic = [1 2 1 2 nan nan nan]; 

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

factors = {'CP', 'Side', 'DC'}; % oly need to label those that 
levels = {{'10' '30' '50'}, {'L' 'R'}, {'1' '3'}}; % oly need to label those that 
levels_aff = {'aff' 'unaff'}; analyse_aff=1;
factor_matrix = [
          1 1 1
          1 1 2
          1 2 1
          1 2 2
          2 1 1
          2 1 2
          2 2 1
          2 2 2
          3 1 1
          3 1 2
          3 2 1
          3 2 2
          ];

if analyse_aff
    levels{2} = levels_aff;
    opt_name = '_aff';
else
    opt_name = '';
end      
      
% create output cell array
header={'Subject'};
for fa = 1:length(factor_matrix)
    header = horzcat(header,[factors{1} levels{1}{factor_matrix(fa,1)} '_' factors{2} levels{2}{factor_matrix(fa,2)} '_' factors{3} levels{3}{factor_matrix(fa,3)}]);
end
acc_all=header;
rt_all=header;
dp_all=header;
cri_all=header;
hit_all=header;
fa_all=header;
dp_all(1,end+1)= {'total'};
cri_all(1,end+1)= {'total'};
hit_all(1,end+1)= {'total'};
fa_all(1,end+1)= {'total'};      

% get data from participant file
subhead = 'Subject';
grphead = 'Group';
inchead = 'Include';
CRPSsidehead = 'CRPSlr';
include_codes = [1];
[~,~,pdata] = xlsread(S.path.datfile);
grp_col = find(strcmp(pdata(1,:),grphead));
sub_col = find(strcmp(pdata(1,:),subhead));
inc_col = find(strcmp(pdata(1,:),inchead));
side_col = find(strcmp(pdata(1,:),CRPSsidehead));
inc_idx = cellfun(@(x) ismember(x,include_codes), pdata(2:end,inc_col), 'UniformOutput', 0);
inc_idx = find(cell2mat(inc_idx));
CRPSsides = pdata(2:end,side_col);
%CRPSsides=CRPSsides(inc_idx);

% add random sides to the healthy subjects (who have no CRPS side)
sidenan = cell2mat(cellfun(@(x) any(x), cellfun(@(x) isnan(x), CRPSsides, 'UniformOutput', 0), 'UniformOutput', 0));
sidepool = CRPSsides(~sidenan);
sn_idx=find(sidenan);
for s=1:length(sn_idx)
    CRPSsides(sn_idx(s)) = sidepool(randi([1 numel(sidepool)]));
end


for d = 1:length(D)
    
    %for op = 1%:length(D(d).Output)
        
        design=D(d).dt.design;
        u = design(3,:); % change = 1, no change = 0;
        for i = 2:size(design,2) % ignore first trial as unlikely to incur a response
            if design(2,i)==0 % mark changes of block as being a stimulus change
                u(1,i)=1;
            end
        end
        
        %comb = [u'];

%         if S.RT.on
%             comb = [comb;D(d).RT.Output(1,:)];
%         end

        % swap sides for people who are right affected
        [hand,dc,cp,bi,blockii,btypes] = blocktype(D(d).dt);
        if analyse_aff && strcmp(CRPSsides{S.subj_pdat_idx(d)},'R')
            i1 = hand==1;
            i2 = hand==2;
            hand(i1)=2;
            hand(i2)=1;
            temp = cell2mat(btypes(:,2));
            i1 = temp==1;
            i2 = temp==2;
            temp(i1)=2;
            temp(i2)=1;
            btypes(:,2) = num2cell(temp)';
            i1 = bi==1;
            i2 = bi==2;
            i3 = bi==3;
            i4 = bi==4;
            i5 = bi==5;
            i6 = bi==6;
            i7 = bi==7;
            i8 = bi==8;
            i9 = bi==9;
            i10 = bi==10;
            i11 = bi==11;
            i12 = bi==12;
            bi(i1)=3;
            bi(i2)=4;
            bi(i3)=1;
            bi(i4)=2;
            bi(i5)=7;
            bi(i6)=8;
            bi(i7)=5;
            bi(i8)=6;
            bi(i9)=11;
            bi(i10)=12;
            bi(i11)=9;
            bi(i12)=10;
        end
        
        % U2: INPUT TYPES
        % create labels for block transitions between hands and within hands
        u2=zeros(1,size(u,2));
        for i = 2:length(blockii) 
            if hand(i) ~= hand(i-1) % label transitions between hands as 7
                u2(blockii(i))=ic(7); 
            elseif dc(i) ~= dc(i-1) % label transitions within hands (from DC1 to DC3) as 5 for left, 6 for right
                if hand(i)==1
                    u2(blockii(i))=ic(5); 
                elseif hand(i)==2
                    u2(blockii(i))=ic(6);
                end
            end
        end
        u(1,blockii)=u2(1,blockii);
        for i = 1:length(u2)
            if design(2,i)~=0
                ii = find(any(cell2mat(btypes(:,1))==design(2,i),2));
                handval = btypes{ii,2};
                dcval = btypes{ii,3};
                if handval==1 && dcval==1
                    u2(i)=ic(1);
                elseif handval==1 && dcval==2
                    u2(i)=ic(2);
                elseif handval==2 && dcval==1
                    u2(i)=ic(3);
                elseif handval==2 && dcval==2
                    u2(i)=ic(4);
                end
            end
        end

        % CHANGE REMAINING ZEROS IN U2 to be same as previous trial
        zi = find(u2==0);
        zi(zi==1)=[]; % ignore first trial
        u2(zi) = u2(zi-1);
        u2(1) = u2(2);

        % binarise inputs
        u(1,:) = u(1,:)>0;
        u(2,:) = u2;
        u(1,isnan(u(2,:)))=NaN;
        D(d).HGF(1).u=u';

        transi = ismember(u2,[1:4]); %use only transitions of interest (1-4 from u2)
        % create index of cond numbers corresponding to each trial
        condi =nan(1,length(u2));
        blockii_end = [blockii length(u2)+1];
        for b = 1:length(bi)
            condi(blockii_end(b):blockii_end(b+1)-1) = bi(b);
        end
        conds = sort(unique(bi));

        % RESPONSES
        if S.RT.on
            if S.fitsim==1
                resp=D(d).RT.Output(1,:);
            elseif S.fitsim==2
                resp = D(d).Output; % for CORE this is always simulated
            end

            if S.fitsim==1
                ISI = 1.0;
                thresh = S.RT.min; % min RT that is realistic
                num_trial_lag = [1 1 1 1]; % number of trials over which response is allowed to lag for each of 0.1 prob, 0.3 prob, 0.5 prob.

                for i = 1:size(design,2)-1
                    if any(cond1 == design(2,i)) || any(cond2 == design(2,i)) || any(cond3 == design(2,i)) || any(cond4 == design(2,i)) 
                        lag = num_trial_lag(1);
                    elseif any(cond5 == design(2,i)) || any(cond6 == design(2,i)) || any(cond7 == design(2,i)) || any(cond8 == design(2,i)) 
                        lag = num_trial_lag(2);
                    elseif any(cond9 == design(2,i)) || any(cond10 == design(2,i)) || any(cond11 == design(2,i)) || any(cond12 == design(2,i)) 
                        lag = num_trial_lag(3);
                    elseif design(2,i)==0;
                        lag = num_trial_lag(4);
                    end

                    resp(2,i)=0; % row 3 is the corrected response time on the corrected trial
                    if u(1,i)>0 % if the stimulus actually changed, but the response occured on the next trial, count it as occuring on the current trial but with a longer RT
                        for j = 1:lag+1
                            move=0;
                            if resp(2,i)==0 && resp(1,i+(j-1))>0 && resp(2,i-1)~=ISI*1+resp(1,i+(j-1)) && resp(1,i+(j-1))+(j-1)*ISI>thresh
                                if u(1,i+(j-1))==0 || j==1 % if we are considering the current trial only, or if there is no stimulus on the next trial, then it's always ok to register the response on te current trial.
                                    move=1;
                                elseif u(1,i+(j-1))>0 && j>1 %if there is a stimulus on the next trial AND
                                    if resp(1,i+(j-1))>0 && resp(1,i+(j-1))<thresh % the response is less than threshold OR
                                        move=1;
                                    elseif (u(1,i+j)==0 && resp(1,i+j)>0) % there is no stim on the subsequent trial but a response
                                        move=1;
                                    end
                                end
                            end
                            if move==1
                                resp(2,i)=resp(1,i+(j-1)) + (j-1)*ISI; % adds on the ISI from RT on next trial and places it on the current trial
                            end
                        end
                    end
                    if u(1,i)==0 && resp(1,i)>0 && i>1 % if no stimulus change, but a response, 
                        if ~any(u(1,i-(j-1):i)>0) % and there was no stim change on previous trial,, so it can't be a delayed response to a stim change
                            if resp(1,i)>thresh;
                                resp(2,i)=resp(1,i); % count it as a false positive response
                            elseif resp(2,i-1)==0
                                resp(2,i-1)=resp(1,i)+1*ISI; % or a response on the previous trial if too fast for current trial
                            end
                        end
                    end
                end 
            end
            if S.fitsim==1
                y = double(resp(2,:)>0); % BINARY response
                y(2,:) = resp(2,:); % RTs
                rtnan = y(2,:);
                rtnan(rtnan==0)=nan;
                y(2,:) = rtnan;
                y(2,:)=log(y(2,:)); % log RTs
                D(d).HGF.y=y';
            elseif S.fitsim==2
                for op = 1:length(resp)
                    y = double(resp(op).pressbutton>0)';
                    D(d).HGF(op).y=y';
                end
                resp=y;
                resp(2,:)=0;
            end

            for op = 1:length(D(d).HGF)
                pc_corr = nan(length(conds),1);
                rt_corr = nan(length(conds),1);
                for c = conds
                    uyr=[];
                    uyr(1,:) = u(1,condi==c);
                    uyr(2,:) = y(1,condi==c);
                    uyr(3,:) = resp(2,condi==c);
                    uyr(4,:) = transi(condi==c);
                    uyr = uyr(:,uyr(4,:)==1);

                    pc_corr(c) = 100*sum(double(uyr(1,:)==uyr(2,:)))/size(uyr,2);
                    hit(c) = length(intersect(find(uyr(1,:)==1),find(uyr(2,:)==1)))/length(find(uyr(1,:)==1)); % hit rate for SDT
                    fa(c) = length(intersect(find(uyr(1,:)==0),find(uyr(2,:)==1)))/length(find(uyr(1,:)==0)); % false alarm rate for SDT
                    rt_corr(c) = mean(uyr(3,intersect(find(uyr(1,:)==1), find(uyr(2,:)==1))));
                    %hit(c+1) = length(intersect(find(u(:,1)==1),find(y==1)))/length(find(u(:,1)==1)); % hit rate for SDT
                    %fa(c+1) = length(intersect(find(u(:,1)==0),find(y==1)))/length(find(u(:,1)==0)); % false alarm rate for SDT
                end
                [dp,cri] = signal_detection(hit,fa);
                D(d).Processed(op).pc_corr = pc_corr;
                D(d).Processed(op).rt_corr = rt_corr;
                D(d).Processed(op).hit = hit;
                D(d).Processed(op).fa = fa;
                D(d).Processed(op).dp = dp;
                D(d).Processed(op).cri = cri;
            end

            if S.save.tables
                acc_all(1+d,1)= pdata(1+S.subj_pdat_idx(d),sub_col);
                rt_all(1+d,1)= pdata(1+S.subj_pdat_idx(d),sub_col);
                dp_all(1+d,1)= pdata(1+S.subj_pdat_idx(d),sub_col);
                cri_all(1+d,1)= pdata(1+S.subj_pdat_idx(d),sub_col);
                hit_all(1+d,1)= pdata(1+S.subj_pdat_idx(d),sub_col);
                fa_all(1+d,1)= pdata(1+S.subj_pdat_idx(d),sub_col);

                acc_all(d+1,2:size(pc_corr,1)+1)=num2cell(pc_corr');
                rt_all(d+1,2:size(rt_corr,1)+1)=num2cell(rt_corr');
                dp_all(d+1,2:size(dp,2)+1)=num2cell(dp');
                cri_all(d+1,2:size(cri,2)+1)=num2cell(cri');
                hit_all(d+1,2:size(dp,2)+1)=num2cell(hit');
                fa_all(d+1,2:size(cri,2)+1)=num2cell(fa');
            end

%     elseif strcmp(part,'part4')
%         load(fullfile(cosdir,cosdirext{1},['CORE' C{2} cosname]));
%         y=nan(length(u),3);
%         
%         %mm=class_proj(1,:)';
%         
%         % standardise mm
%         %mm=zscore(mm);
%         
%         % add raw mm to column 3 of y
%         try
%             y(tnums,3)=mm;
%         catch
%             y(tnums,3)=MM.mm;
%         end
%         if mm_trials
%             y(:,3) = y(:,3).*u(:,1);
%         end
%         if mm_positive
%             y(:,3) = y(:,3).*(y(:,3)>0);
%         end
%         y(y(:,3)==0,3) = nan;
%         
%         % class predictions for softmax function
%         %mm_thresh = 0; % number of SDs
%         %class1=1*single(mm>mm_thresh);%CAB 
%         %class2=0*single(mm<=mm_thresh);%CAB 
%         %predicted=class1+class2;%CAB 
%         %y(tnums,4)=predicted;

        %D(d).Processed.comb = comb;
        %end
    %end
    
    
        % average simualted data
        if S.meansim && S.fitsim==2
            S.mean_fields = {};
            if S.accuracy.on
                S.mean_fields = [S.mean_fields {'pc_corr','rt_corr','hit','fa','dp','cri'}];
            end
            for mf = 1:length(S.mean_fields)
                if iscell(D(d).Processed(1).(S.mean_fields{mf}))
                    for i1 = 1:length(D(d).Processed(1).(S.mean_fields{mf}))
                        if iscell(D(d).Processed(1).(S.mean_fields{mf}){i1})
                            for i2 = 1:length(D(d).Processed(1).(S.mean_fields{mf}){i1})
                                temp_mf=[];
                                for op = 1:length(D(d).Processed)
                                    temp_mf(op,:) = D(d).Processed(op).(S.mean_fields{mf}){i1}{i2};
                                end
                                D(d).Processed(1).(S.mean_fields{mf}){i1}{i2} = mean(temp_mf,1);
                            end
                        else
                            temp_mf=[];
                            for op = 1:length(D(d).Processed)
                                temp_mf(op,:) = D(d).Processed(op).(S.mean_fields{mf}){i1};
                            end
                            D(d).Processed(1).(S.mean_fields{mf}){i1} = mean(temp_mf,1);
                        end
                    end
                else
                    temp_mf=[];
                    for op = 1:length(D(d).Processed)
                        temp_mf(op,:) = D(d).Processed(op).(S.mean_fields{mf});
                    end
                    D(d).Processed(1).(S.mean_fields{mf}) = mean(temp_mf,1);
                end
            end
            D(d).Processed = D(d).Processed(1);
        end
    
    end
    
end


if S.save.tables
    if S.fitsim==1
        spath=S.path.prep;
    elseif S.fitsim==2
        spath=fullfile(S.path.prep,'sim');
    end
    cd(spath)
    sname = datestr(now,30);
    xlswrite(['accuracy' opt_name '_' sname '.xlsx'],acc_all);
    xlswrite(['reactiontime' opt_name '_' sname '.xlsx'],rt_all);
    xlswrite(['dprime' opt_name '_' sname '.xlsx'],dp_all);
    xlswrite(['criterion' opt_name '_' sname '.xlsx'],cri_all);
    xlswrite(['hitrate' opt_name '_' sname '.xlsx'],hit_all);
    xlswrite(['farate' opt_name '_' sname '.xlsx'],fa_all);
end


