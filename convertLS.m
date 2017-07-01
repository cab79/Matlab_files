function convertLS(D)

% CONVRSION FROM LONG TO SHORT FORMAT AND VICE VERSA
% "Long" format data has columns for different factors in the leftmost columns,
% and a single rightmost column containing the data.
% "Short" format data can have some or no factor columns, and has a number of rightmost columns
% containing data that are combinations of factors. 

% To convert between 2 different "Short" formats, the script needs to be run twice:
% 1. Short to Long
% 2. Long to new Short

% load file
[~,~,pdata] = xlsread(D.pdatfile);

% sort
if strcmp(D.format,'long') && ~isempty(D.fact_rm) % must sort by factors other than those being removed FIRST
    truesortin = D.sortby_in;
    D.sortby_in = 1:D.Nfact_in;
    D.sortby_in(D.fact_rm)=[];
    if intersect(truesortin,D.fact_rm)
        D.sortby_in = [D.sortby_in intersect(truesortin,D.fact_rm)]; % must sort by factors being removed LAST
    end
end
for s = 1:length(D.sortby_in)
    fh = @(x) all(isnan(x(:)));
    nanind=(cellfun(fh, pdata(:,s)));
    try pdata(nanind,s) = {''};
    catch
    end
end
if ~isempty(D.sortby_in); pdata(2:end,:) = sortrows(pdata(D.Nhead+1:end,:),D.sortby_in);end

% collect data
dat_in = cell2mat(pdata(D.Nhead+1:end,D.Nfact_in+1:D.Nfact_in+D.Ndat_in));
dat_head_in = pdata(D.Nhead,D.Nfact_in+1:D.Nfact_in+D.Ndat_in);
dat_head_in = cellfun(@num2str,dat_head_in,'UniformOutput',0);
Sdat_in = size(dat_in);

% collect factor info
fact_in = pdata(1:D.Nhead,1:D.Nfact_in);
fact_val_in = {};
fact_uniq_in = {};
Nlev = [];
for f = 1:length(fact_in)
    if isnumeric(pdata{D.Nhead+1,f})
        fact_val_in{f} = cell2mat(pdata(D.Nhead+1:end,f));
    else
        fact_val_in{f} = pdata(D.Nhead+1:end,f);
        fh = @(x) all(isnan(x(:)));
        try fact_val_in{f}(cellfun(fh, fact_val_in{f})) = {''};
        catch
        end
    end
   
    fact_uniq_in{f} = unique(fact_val_in{f});
    Nlev_in(f) = length(fact_uniq_in{f});
end

% process outputs
pdataout = {};
fact_val_out ={};
if strcmp(D.format,'short') % short format; convert to long
    
    % reshape input data to one column output
    dat_out = reshape(dat_in',Sdat_in(1)*Sdat_in(2),1);
    
    % list of output factors
    fact_out = horzcat(fact_in,D.fact_new);
    Nfact_out = length(fact_out);
    
    % values of output factors
    for f = 1:length(fact_in)
        fact_val_out{f} = reshape(repmat(fact_val_in{f},1,Sdat_in(2))',Sdat_in(1)*Sdat_in(2),1);
    end
    
    % values of new output factors
    fact_val_new = {};
    for f = 1:length(D.fact_new)
        Nrep1 = D.Nlev_fact_new(f:end);
        Nrep1(1)=[];
        Nrep2 = D.Nlev_fact_new(1:f);
        Nrep2(end)=[];
        if length(dat_head_in)==D.Nlev_fact_new
            lev_suff = dat_head_in; 
        else
            lev_suff = cellfun(@num2str,num2cell(1:D.Nlev_fact_new),'UniformOutput',0);
        end
        for fl = 1:D.Nlev_fact_new(f)
            fact_val_new{f}{fl,1} = [D.fact_new{f} lev_suff{fl}];
        end
        if ~isempty(Nrep1); fact_val_new{f} = reshape(repmat(fact_val_new{f},1,prod(Nrep1))',Sdat_in(2)/prod(Nrep2),1);end
        if ~isempty(Nrep2); fact_val_new{f} = reshape(repmat(fact_val_new{f},prod(Nrep2),1)',Sdat_in(2),1);end
        fact_val_new{f} = repmat(fact_val_new{f},Sdat_in(1),1);
    end
    fact_val_out = horzcat(fact_val_out,fact_val_new);
    
    Ndat_head = {'Data'};
    
    savename = 'long';
elseif strcmp(D.format,'long') % long format; convert to short
    
    % list of output factors and values
    fact_out = fact_in;
    fact_out(D.fact_rm)=[];% remove factors
    Nfact_out = length(fact_out);
    fact_val_out = fact_val_in; % create factor values from input value
    fact_val_out(D.fact_rm)=[];% remove factor values
    
    % number of combinations of factors being removed = number of output data columns
    Ndat_out=1;
    for rf = 1:length(D.fact_rm)
        Ndat_out = Ndat_out*Nlev_in(D.fact_rm(rf));
    end
    
    % reshape input data to multiple column output
    dat_out = reshape(dat_in,Ndat_out,Sdat_in(1)/Ndat_out)';

    % reshape output factors
    for f = 1:length(fact_val_out)
        if isnumeric(fact_val_out{f}(1))
            fact_val_out{f} = num2cell(reshape(cell2mat(fact_val_out(f)),Ndat_out,Sdat_in(1)/Ndat_out)');
        else
            fact_val_out{f} = reshape(fact_val_out{f},Ndat_out,Sdat_in(1)/Ndat_out)';
        end
        fact_val_out{f} = fact_val_out{f}(:,1);
    end
    
    fact_val_rm = pdata(D.Nhead+1:D.Nhead+Ndat_out,D.fact_rm);
    Ndat_head = {};
    for nc = 1:Ndat_out
        for f = 1:length(D.fact_rm)
            if length(Ndat_head) < nc; Ndat_head{nc,1} = '';end
            if isnumeric(fact_val_rm{nc,f}); 
                dat_head_name = num2str(fact_val_rm{nc,f});
            else
                dat_head_name = fact_val_rm{nc,f};
            end
            Ndat_head{nc,1} = [Ndat_head{nc,1} '_' fact_in{D.fact_rm(f)} '_' dat_head_name];
        end
        Ndat_head{nc,1} = Ndat_head{nc,1}(2:end);
    end
    
    savename = 'short';
end

% build output file
Ndat_out = size(dat_out);
pdatout(1,1:length(fact_out)) = fact_out;
pdatout(1,length(fact_out)+1:length(fact_out)+length(Ndat_head)) = Ndat_head;
for f = 1:length(fact_out)
    if isnumeric(fact_val_out{f}(1))
        pdatout(D.Nhead+1:D.Nhead+Ndat_out(1),f) = num2cell(fact_val_out{f});
    else
        pdatout(D.Nhead+1:D.Nhead+Ndat_out(1),f) = fact_val_out{f};
    end
end
pdatout(D.Nhead+1:D.Nhead+Ndat_out(1),f+1:f+Ndat_out(2)) = num2cell(dat_out);

% sort
if ~isempty(D.sortby_out); pdatout(D.Nhead+1:end,:) = sortrows(pdatout(D.Nhead+1:end,:),D.sortby_out);end;

if ~isempty(D.avg)
    newdatout_temp = pdatout;
    %if isnumeric(pdatout{D.Nhead+1,D.avg})
    %    avg_col = cell2mat(pdatout(D.Nhead+1:end,D.avg));
    %    uniqlev = unique(avg_col);
    %    rowind={};
    %    for ul = 1:length(uniqlev)
    %        rowind{ul,1} = find(avg_col==uniqlev(ul));
    %    end
    %else
    avgind = [];
    for na = 1:length(D.avg)
        if isnumeric(pdatout{D.Nhead+1,D.avg(na)})
            [~,~,avgind(:,na)] = unique(cellfun(@num2str,pdatout(D.Nhead+1:end,D.avg(na)),'uni',0),'stable'); %unique indices
        else
            [~,~,avgind(:,na)] = unique(pdatout(D.Nhead+1:end,D.avg(na)),'stable'); %unique indices
        end
    end
    [rowuniq,firstrowuniq] = unique(avgind,'rows','stable');%unique row indices 
    uniqlev = pdatout(D.Nhead+firstrowuniq,D.avg);
    rowind=cell(1,size(uniqlev,1));
    for ur = 1:size(uniqlev,1) % for each row of unique levels
        rowcolind={};
        for c = 1:size(uniqlev,2) % for each column of unique levels
            % find rows of fact_val_out
            if isnumeric(uniqlev{ur,c})
                rowcolind{c} = find(cell2mat(pdatout(D.Nhead+1:end,D.avg(c)))==uniqlev{ur,c});
            else
                rowcolind{c} = find(strcmp(pdatout(D.Nhead+1:end,D.avg(c)),uniqlev{ur,c}));
            end
            if isempty(rowind{ur}); rowind{ur} = rowcolind{c}; end;
            rowind{ur}=intersect(rowind{ur},rowcolind{c});
        end
    end
    
    rmfact = 1:Nfact_out;
    rmfact(D.avg) = [];
    newdatout_temp(:,rmfact) = [];
    newdatout = newdatout_temp(1,:);
    for nc = 1:size(newdatout_temp,2)
        if isnumeric(newdatout_temp{D.Nhead+1,nc})
            meandat=[];
            for ul = 1:length(uniqlev)
                dat = cell2mat(newdatout_temp(D.Nhead+1:end,nc));
                meandat(ul,1) = nanmean(dat(rowind{ul}));
            end
            newdatout(D.Nhead+1:D.Nhead+length(meandat),nc) = num2cell(meandat);
        else
            meandat={};
            for ul = 1:length(uniqlev)
                dat = newdatout_temp(D.Nhead+1:end,nc);
                meandat{ul,1} = dat{rowind{ul}(1)};
            end
            newdatout(D.Nhead+1:D.Nhead+length(meandat),nc) = meandat;
        end
        
    end
    pdatout = newdatout;
end

[pth,nme,ext] = fileparts(D.pdatfile);
fname = fullfile(pth,[nme '_' savename '_' datestr(datetime,30) ext]);
xlswrite(fname,pdatout);