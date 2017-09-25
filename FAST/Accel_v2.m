%% create settings structure (don't change)
clear all
S=struct;
dbstop if error

%% INPUT FILE SETUP

% accelerometer data directory
S.acc_dir = 'C:\Matlab\FAST';
% accelerometer data generic file name (without subject identifier)
S.acc_file = '_WTV_MS_1.xlsx';
% column header for date and time
S.acc_datehead = {'date','epoch'};
% column header for data to be analysed
S.acc_data = 'steps';

% pain ratings data directory
S.pr_dir = 'C:\Matlab\FAST';
% pain ratings data generic file name (without subject identifier)
S.pr_file = '_PainScores_MS.xlsx';
% column header for date and time
S.pr_datehead = 'Date and time';
% column header for pain rating
S.pr_data = 'Pain score';

%% OUTPUT FILE DIRECTORY AND NAME
S.out_dir = 'C:\Matlab\FAST';
% save new file for each subject (saveeach = 1), or all subjects in one
% file (=0)? (0 is not currently programmed)
S.saveeach = 1;
% generic file name (subject identifier will be added as prefix)
S.out_file = '_output.xlsx'; 

%% ANALYSIS SETTINGS

% First, define the time windows over which to apply functions (e.g. functions are simple summary measures such as the mean) to the data
% Also, define the accelerometer time window with respect to the pain ratings time window (S.acc_cen)

% time window to apply to accelerometer data
S.acc_tw = '2'; % current options: 'daily'. Other options: a number of hours, e.g. '2' 
% time window to apply to pain ratings data
S.pr_tw = '1'; % current options: 'daily'. Other options: a number of hours, e.g. '1' 
% centering of accelerometer data with respect to pain ratings outputs.
% if left empty as '', analysis of accelerometer data will proceed
% independently of pain ratings data, or will overlap in the case of
% 'daily' options above. Otherwise, a number e.g.:
% '-2' = S.acc_tw starts two hours BEFORE the pain rating
% '2' = S.acc_tw starts two hours AFTER the pain rating
S.acc_cen = '-2'; % only used for options above other than 'daily'. 

% update output name
S.out_file = ['_accTW_' S.acc_tw '_prTW_' S.pr_tw '_accCEN_' S.acc_cen S.out_file];

% List functions to apply to each type of data
% Output file will include the results of these functions

% ACCEL FUNCTION 1
% name with no spaces
S.acc_fun(1).name = 'mean'; 
% actual function using Matlab syntax. X = N-dimensional data array.
S.acc_fun(1).fun = 'mean(X)'; 

% ACCEL FUNCTION 2
% name with no spaces
S.acc_fun(2).name = 'median'; 
% actual function using Matlab syntax. X = N-dimensional data array.
S.acc_fun(2).fun = 'median(X)'; 

% ACCEL FUNCTION 3
% name with no spaces
S.acc_fun(3).name = 'std'; 
% actual function using Matlab syntax. X = N-dimensional data array.
S.acc_fun(3).fun = 'std(X)'; 

% ACCEL FUNCTION 4
% name with no spaces
S.acc_fun(4).name = 'range'; 
% actual function using Matlab syntax. X = N-dimensional data array.
S.acc_fun(4).fun = 'max(X)-min(X)'; 

% ACCEL FUNCTION 5
% name with no spaces
S.acc_fun(5).name = 'sum'; 
% actual function using Matlab syntax. X = N-dimensional data array.
S.acc_fun(5).fun = 'sum(X)'; 

% PAIN RATINGS FUNCTION 1
% name with no spaces
S.pr_fun(1).name = 'mean'; 
% actual function using Matlab syntax. X = N-dimensional data array.
S.pr_fun(1).fun = 'mean(X)'; 

%% RUN SCRIPT

% create date/time string to add to output file name (not used)
%runtime = datestr(now,30);
% list all accelerometer files in the directory
afiles = dir(fullfile(S.acc_dir,['*' S.acc_file]));

% create output headers defining dates/times
D = {'Date','Time'};
nHead = size(D);

for f = 1:length(afiles)
    
    clear pr_head pr_out acc_head acc_out
    
    % load each accel file in turn and also load it's corresponding pain
    % ratings data file (if there is one)
    
    % obtain file name
    acc_file = afiles(f).name;
    % identify subject
    temp = strsplit(acc_file,'_');
    subname = [temp{1} '_' temp{2}];
    % obtain corresponding pain ratings file name
    pr_file = [subname S.pr_file];
    % load each data file
    [~,~,acc] = xlsread(fullfile(S.acc_dir,acc_file));
    try
        [~,~,pr] = xlsread(fullfile(S.pr_dir,pr_file));
    catch
        % if there are no pain ratings: create empty variables
        pr=[];
        S.acc_cen = ''; % cannot centre acc data with respect to pr
    end
    % housekeeping: remove first header
    acc(1,:)=[];
    
    %% get dates/times
    
    % for accel data:
    % column of file with date/time
    acc_dt_col = ismember(acc(1,:),S.acc_datehead);
    % datetime data
    acc_dt_dat = acc(2:end,acc_dt_col);
    % parse into dates and times
    acc_ds = acc_dt_dat(:,1);
    acc_ds = cellstr(datestr(datenum(acc_ds,'dd/mm/yyyy'),29));
    acc_ts = acc_dt_dat(:,2);
    acc_ts = cellfun(@(x) datestr(x,13),acc_ts,'UniformOutput',0);
    % create datestr
    acc_dt_str = strcat(acc_ds, {' '}, acc_ts);
    % create datevec
    acc_dt_vec = datevec(acc_dt_str,'yyyy-mm-dd HH:MM:SS');
    
    % for PR data:
    if ~isempty(pr)
        % column of file with date/time
        pr_dt_col = ismember(pr(1,:),S.pr_datehead);
        % datetime data
        pr_dt_dat = pr(2:end,pr_dt_col);
        % fix Excel bug that causes 00:00 times to not be imported via xlsread
        pr_dt_dat(cellfun(@length,pr_dt_dat)~=19)=strcat(pr_dt_dat(cellfun(@length,pr_dt_dat)~=19),' 00:00:00');
        % find indices of first unique entry and remove duplicates
        [~,pr_dt_uni,~] = unique(pr_dt_dat,'stable');
        pr_dt_dat = pr_dt_dat(pr_dt_uni);
        % convert to datevec
        pr_dt_vec = datevec(pr_dt_dat,'dd/mm/yyyy HH:MM:SS');
        % parse into dates and times
        pr_ds = datestr(pr_dt_vec,29);
        pr_ts = datestr(pr_dt_vec,13);
        %temp = cellfun(@(x) strsplit(x,' '),pr_dt_dat,'UniformOutput',0);
        %pr_ds = cellfun(@(x) x(1),temp,'UniformOutput',0);
        %pr_ts = cellfun(@(x) x(2),temp,'UniformOutput',0);
        %pr_ds=[pr_ds{:}]'; % dates
        %pr_ds = cellstr(datestr(datenum(pr_ds,'dd/mm/yyyy'),29));
        %pr_ts=[pr_ts{:}]'; % times
    end
    
    
    %% analyse pain ratings data (if the file exists)
    if ~isempty(pr)
        % column of file with pain ratings data
        pr_col = strmatch(S.pr_data,pr(1,:));
        % pain ratings data
        pr_dat = pr(2:end,pr_col);
        % values only from unique dates/times
        pr_dat = pr_dat(pr_dt_uni);
        
        % find indices of data according to the specified time window
        if strcmp(S.pr_tw,'daily')
            % unique dates and their index
            [pr_tu,~,pr_ti] = unique(pr_ds,'stable');
            % time column is empty
            pr_tu(:,2)=cell(length(pr_tu),1);
            pr_tu(:,2)={''};
        elseif ~isnan(str2double(S.pr_tw)) % if a scalar value, e.g. number of hours
            % run through times and chunk together if within S.pr_tw of
            % each other
            pr_ti = zeros(length(pr_dt_vec),1); % assignment of times to a chunk index
            for i = 1:length(pr_dt_vec)-1
                % if the difference between two adjacent times is less than
                % S.pr_tw
                if abs(etime(pr_dt_vec(i+1,:),pr_dt_vec(i,:)))/60 < str2double(S.pr_tw)*60
                    try 
                        % if the current time is already part of a chunk
                        if pr_ti(i)==pr_ti(i-1)
                            % assign the next time also to the same chunk
                            pr_ti(i+1) = pr_ti(i);
                        else
                            % otherwise create a new chunk
                            pr_ti(i:i+1) = max(pr_ti)+1;
                        end
                    catch
                        % by default, create a new chunk
                        pr_ti(i:i+1) = max(pr_ti)+1;
                    end
                elseif pr_ti(i)==0
                    pr_ti(i) = max(pr_ti)+1;
                end
            end
            % complete the last value if it's not part of a chunk already
            if pr_ti(end)==0
                pr_ti(end) = max(pr_ti)+1;
            end
            
            % create new times that are the centre (mean) of each chunk
            for i = unique(pr_ti)'
                pr_tu_vec(i,:) = datevec(mean(datenum(pr_dt_vec(pr_ti==i,:))));
            end
            
            % turn pr_tu into 2-column cellstr array of date/time
            pr_ds = datestr(pr_tu_vec,29);
            pr_ts = datestr(pr_tu_vec,13);
            pr_tu = cellstr(pr_ds);
            pr_tu(:,2) = cellstr(pr_ts);
        elseif strcmp(S.pr_tw,'all')
            pr_tu = {'all',''};
            pr_ti = ones(length(pr_ds),1);
        else
            error('this time window is not yet programmed!');
        end
        
        % calculate functions over specified time window
        for fi = 1:length(S.pr_fun)
            % create column header
            pr_head(1,fi) = {[genvarname(S.pr_data) '_' genvarname(S.pr_fun(fi).name)]};
            % values of pr_ti to use as indices on pr_dat
            val = unique(pr_ti)';
            val(val==0)=[];
            for i = val
                % get data
                X = cell2mat(pr_dat(pr_ti==i,:));
                % apply function
                eval(['fout = ' S.pr_fun(fi).fun ';']);
                pr_out(val==i,fi) = fout;
            end
        end 
        
    else
        pr_head = {};
        pr_out = [];
    end
    
    %% analyse accel data
    
    % column of file with accel data
    acc_col = find(ismember(acc(1,:),S.acc_data));
    % accel data
    acc_dat = acc(2:end,acc_col);

    % find indices of data according to the specified time window
    if strcmp(S.acc_tw,'daily')
        % unique dates and their index
        [acc_tu,~,acc_ti] = unique(acc_ds,'stable');
        % time column is empty
        acc_tu(:,2)=cell(length(acc_tu),1);
        acc_tu(:,2)={''};
    elseif ~isnan(str2double(S.acc_tw)) % if a scalar value, e.g. number of hours
        % if the data is to be analysed according to timings of pain ratings
        if ~isempty(pr)
            % find time windows (start and end times) defining S.acc_tw,
            % centred on S.acc_cen
            acc_tw = [];
            for i = 1:length(pr_tu) % unique PR times
                s_num = addtodate(datenum(pr_tu_vec(i,:)), str2double(S.acc_cen), 'hour'); % start vector
                e_num = addtodate(datenum(pr_tu_vec(i,:)), str2double(S.acc_cen)+str2double(S.acc_tw), 'hour'); % end vector
                acc_tw(i,:) = [s_num,e_num];
            end
            acc_tu=pr_tu; % only if there is some acc data in that tw?
            
        % otherwise, if there are no timings available from pain ratings
        else
            % identify first and last datetime
            first_num = datenum(acc_dt_vec(1,:));
            last_num = datenum(acc_dt_vec(end,:));
            % create list of datetimes between start and end date chunked
            % by S.acc_tw
            e_num = 0;
            ci = 0;
            acc_tw = [];
            while e_num<last_num
                ci = ci+1;
                s_num = addtodate(first_num, (ci-1)*str2double(S.acc_tw), 'hour'); % start vector
                e_num = addtodate(first_num, ci*str2double(S.acc_tw), 'hour'); % end vector
                acc_tw(ci,:) = [s_num,e_num];
            end
            acc_tu = cellstr(datestr(acc_tw(:,1),29));
            acc_tu(:,2) = cellstr(datestr(acc_tw(:,1),13));
        end
        % run through times and chunk together according to acc_tw
        acc_dt_num = datenum(acc_dt_vec);
        acc_ti = zeros(length(acc_dt_num),1); % assignment of times to a chunk index
        ci=0;
        tu_ind = 1:size(acc_tu,1);
        for i = 1:size(acc_tw,1)
            ind = find(acc_dt_num >= acc_tw(i,1) & acc_dt_num < acc_tw(i,2));
            if isempty(ind)
                % remove dametimes in which there is no data 
                acc_tu(tu_ind==i,:)=[];
                tu_ind(tu_ind==i)=[];
            else
                ci = ci+1;
                acc_ti(ind) = ci;
            end
        end
        % remove any datetimes in acc_tu if they no longer have
        % corresponding chunks of data - can occur if there are two sets of
        % pain ratings that are close together and are used to calculate
        % windows around the exact same accel data - only the second chunk
        % remains in acc_ti
        [incl,~] = ismember(1:length(acc_tu),acc_ti);
        acc_tu = acc_tu(incl,:);
    
    elseif strcmp(S.acc_tw,'all')
        acc_tu = {'all',''};
        acc_ti = ones(length(acc_ds),1);
    else
        error('this time window is not yet programmed!');
    end

    
    % calculate functions over specified time window
    for fi = 1:length(S.acc_fun)
        % create column header
        acc_head(1,fi) = {[genvarname(S.acc_data) '_' genvarname(S.acc_fun(fi).name)]};
        % values of acc_ti to use as indices on acc_dat
        val = unique(acc_ti)';
        val(val==0)=[];
        for i = val
            % get data
            X = cell2mat(acc_dat(acc_ti==i,:));
            % apply function
            eval(['fout = ' S.acc_fun(fi).fun ';']);
            acc_out(val==i,fi) = fout;
        end
    end  
    
    
    %% output dates and data
    if isempty(pr)
        out_tu = acc_tu;
        %acc_tu_i = 1:size(out_tu,1);
        %out_data = acc_out;
        acc_out = num2cell(acc_out);
    else
    % if the PR data exists, this will require interleaving dates/times
    % from both data types.
    % requires matching up the PR and accel dates/times
    % finds unique dates as well as indices of those unique outputs that
    % match the rows of each data type (PR and accel)
        temp = [pr_tu;acc_tu];
        [~,i1,i2] = unique(cell2mat(temp),'rows');
        out_tu = temp(i1,:);
        pr_tu_i = 1:size(pr_tu,1);
        acc_tu_i = size(pr_tu,1)+1:size(pr_tu,1)+size(acc_tu,1);
        pr_tu_i = i2(pr_tu_i);
        acc_tu_i = i2(acc_tu_i);
        
        % output data arrays
        temp = cell(size(out_tu,1),size(acc_out,2));
        temp(acc_tu_i,:) = num2cell(acc_out);
        acc_out = temp;
        
        temp = cell(size(out_tu,1),size(pr_out,2));
        temp(pr_tu_i,:) = num2cell(pr_out);
        pr_out = temp;
    end
    
    %% create output cell array and save
    % add PR data output headers
    allhead = horzcat(pr_head,acc_head);
    D(nHead(1),nHead(2)+1:nHead(2)+length(allhead)) = allhead;
    % add date/time data to output
    D(nHead(1)+1:nHead(1)+size(out_tu,1),1:nHead(2)) = out_tu;
    % add PR data to output
    D(nHead(1)+1:nHead(1)+size(out_tu,1),nHead(2)+1:nHead(2)+size(pr_out,2)) = pr_out;
    % add accel data to output
    D(nHead(1)+1:nHead(1)+size(out_tu,1),nHead(2)+size(pr_out,2)+1:nHead(2)+size(pr_out,2)+size(acc_out,2)) = acc_out;
    % save 
    if S.saveeach
        xlswrite(fullfile(S.out_dir,[subname S.out_file]),D); 
    else
        error('saving multiple subjects in one file not yet implemented');
    end
    
end
    