
%% run
function gplot_waveforms(P)
dbstop if error

if ~isfield(P,'timezero')
    P.timezero = [];
end

% change zero time
if ~isempty(P.timezero)
    P.xlines = P.xlines-P.timezero;
end

%% load EEGLAB data using filenames from an SPM file and reformat for plotting topography
if P.plot_topo && ~exist(fullfile(P.eeglab_path,'ERP_DAT.mat'),'file')
    P.st_string='mspm12_';
    P.en_string='\scond';
    P.savename = 'ERP_DAT.mat';
    gplotprepare_eeglabdata_from_spm(P)
elseif P.plot_topo
    load(fullfile(P.eeglab_path,'ERP_DAT.mat'));
end

%% plot cluster waveforms
% load waveform data
load(fullfile(P.spm_path,P.clusdir,P.wfname));

% load SPM design and list factors
load(fullfile(P.spm_path,P.batch));
P.fact = {matlabbatch{1,1}.spm.stats.factorial_design.des.fblock.fac(:).name};

% identify file indices relating to factor of interest and subject
fact_col = [];
for c = 1:length(P.facplot) % using a loop ensures that fact_col is in correct order for condition labels applied later
    fact_col(c) = find(ismember(P.fact,P.facplot(c)));
end
sub_col = find(ismember(P.fact,P.subfactname));
fact_ind = S.Fm(:,1+fact_col);
sub_ind = S.Fm(:,1+sub_col);
if any(strcmp(P.fact,'Time')) && S.timelev
    time_ind = S.Fm(:,1+find(strcmp(P.fact,'Time')));
    fact_ind = fact_ind(time_ind==S.timelev,:);
    sub_ind = sub_ind(time_ind==S.timelev);
end

if ~isempty(P.plotclus)
    % restrict to selected clusters
    clnames = P.plotclus';
else
    % obtain all cluster names
    clnames = fieldnames(S.wf);
end

% load example cluster image to extract time information
Cnii = load_nii(fullfile(P.spm_path,P.clusdir,[clnames{1} '.nii']));
t_res = Cnii.hdr.hist.srow_z(3);
t_off = Cnii.hdr.hist.srow_z(4);
t_dim = Cnii.hdr.dime.dim(4);
t_last = t_off+t_res*(t_dim-1);
itimes = t_off:t_res:t_last;

% load cluster statistics
load(fullfile(P.spm_path,P.clusdir,'cluster_table.mat'));

% colormap of background polygons
%colormap(gray);
%cmap=colormap;
%close all
range=[];
% for each cluster, 
for cl = 1:length(clnames)
    C=strsplit(clnames{cl},'_');
    cllabel = C{1};
    
    % select cluster data
    wf = S.wf.(clnames{cl});
    % select wf if struct (e.g. source data)
    if isstruct(wf)
        itimes = wf.time*1000;
        wf=wf.wf;
    end
    % identify unique rows of the combination of factors
    [~,Frows,WFrows] = unique([fact_ind sub_ind],'rows');
    
    % create new factor indices for unique rows.
    Fi = fact_ind(Frows,:);
    Si = sub_ind(Frows,:);
    
    % average wf data over non-unique rows
    wff = cell(length(unique(WFrows)),1);
    for i = unique(WFrows)'
        if size(wf{1},1)>size(wf{1},2) % timepoints in first dimension
            wff{i}=mean(cell2mat(wf(WFrows==i)'),2);
        else
            wff{i}=mean(cell2mat(wf(WFrows==i)),1)';
        end
    end
    
    % if there is more than one factor, the first factor levels produce a
    % plot each
    P.ptitle={};
    Fi_ind={};
    if length(P.facplot)>1
        P.Nxplots=length(P.cval{1,2});
        % titles
        for pt = 1:P.Nxplots
            P.ptitle{pt} = P.cval{1,1}{P.cval{1,2}(pt)};
            Fi_ind{pt} = find(Fi(:,1)==P.cval{1,2}(pt)); % indices of Fi (and wf) for each plot
        end
    else
        P.Nxplots = 1;
        P.ptitle{1} = '';
        Fi_ind{1} = 1:length(Fi);
    end
    
    % extract cluster statistics
    ct_ind = find(ismember(clustable(:,1),cllabel));
    if ~isempty(P.poly)
        E_val = P.poly;
    else
        %F_val = clustable{ct_ind,9};
        E_val = [clustable{ct_ind,11:12}];
        P_val = [clustable{ct_ind,10}];
    end
    
    % assign x lines and polygons
    if ~isempty(P.timezero)
        P.poly = E_val-P.timezero;
        P.xlinesolid = P_val-P.timezero;
    else
        P.poly = E_val;
        P.xlinesolid = P_val;
    end
    
%end

% for each cluster, 
%for cl = 1:length(clnames)
    
    for p = 1:P.Nxplots
        % plot data trajectories
        % gramm supports 2D inputs for X and Y data (as 2D array or cell of
        % arrays), which is particularly useful for representing repeated
        % trajectories. The grouping data is then given per trajectory, given
        % as a 1xNtraj cellstr.
        
        % construct condition labels (of last or only factor)
        condlev=P.cval{end,1};
        factind=Fi(Fi_ind{p},end);
        fi{p}=ismember(factind,P.cval{end,2});
        factind = factind(fi{p});
        P.cond{p} = condlev(factind);
        
        % set x and y axis data to plot
        P.y=wff(Fi_ind{p});
        P.y=P.y(fi{p});
        if ~isempty(P.xlimits)
            xlim = dsearchn(itimes',P.xlimits')';% x values of the selected segment to plot
            P.x = itimes(xlim(1):xlim(2));
            for n = 1:length(P.y)
                P.y{n} = P.y{n}(xlim(1):xlim(2));
            end
        else
            P.x=itimes;
        end
        
        % change zero time
        if ~isempty(P.timezero)
            P.x = P.x-P.timezero;
        end
        
        
    end
    
    %-------- plot time series -----------%
    P.gtitle = clnames{cl};
    g = gplot_timeseries(P)
    %-------------------------------------%
    
    if P.save_waveforms
        psize = 18;
        g.export('file_name',[cllabel '_WF'],'export_path',fullfile(P.spm_path,P.clusdir),'file_type','svg','width',P.Nxplots*psize,'height',psize,'units','centimeters');
    end
    
    if P.plot_topo
        % requires:
            % DAT (output from plotprepare_)
            % WFrows
            % P.Nxplots
            % Fi_ind
            % fi
            % cond
            % chanlocs
            % P
            % Eval
            % EEGtimes
        
        % average ERPs over non-unique rows
        DATa = cell(length(unique(WFrows)),1);
        for i = unique(WFrows)'
            DATa{i}=mean(cat(3,DAT{WFrows==i}),3);
        end
        f1=figure
        pln=0;
        for p = 1:P.Nxplots
            y=DATa(Fi_ind{p});
            y=y(fi{p});
            [~,~,condind]=unique(P.cond{p},'stable');
            unicond = unique(condind)';
            for cn = unicond
                peakdata=y(condind==cn);
                peakdata = mean(cat(3,peakdata{:}),3);
                plotchans=1:length(chanlocs);
                plotchans(P.no_plot_ele)=[];
                %[~,markchans] = intersect(plotchans,tp);
                if any(P.topo_subtimewin) && length(P.topo_subtimewin)==1 % multiple plots within a range
                    cluswin = E_val(end)-E_val(1);
                    nlat = floor(cluswin/P.topo_subtimewin);
                    if nlat>1
                        lats = linspace(E_val(1),E_val(end),nlat);
                        lats = dsearchn(EEGtimes',lats');
                    else
                        lats = [find(EEGtimes==E_val(1)) find(EEGtimes==E_val(end))];
                    end
                    
                    for ln = 1:length(lats)-1
                        pln = pln+1;
                        lat=lats(ln:ln+1);
                        subplot(P.Nxplots+length(unicond),length(lats)-1,pln); topoplot(mean(peakdata(:,lat(1):lat(2)),2), chanlocs,'maplimits','absmax','electrodes','on','plotchans',plotchans);%,'emarker2',{markchans,'o','w',7,1}); 
                        title(num2str(EEGtimes(lat)'))
                    end
                elseif any(P.topo_subtimewin) && length(P.topo_subtimewin)==2 % specified time window
                    pln = pln+1;
                    lat = dsearchn(EEGtimes',P.topo_subtimewin');
                    subplot(P.Nxplots+length(unicond),P.Nxplots,pln); topoplot(mean(peakdata(:,lat),2), chanlocs,'maplimits','absmax','electrodes','on','plotchans',plotchans);%,'emarker2',{markchans,'o','w',7,1}); 
                    title(num2str(EEGtimes(lat)'))
                else
                    pln = pln+1;
                    lat = find(EEGtimes==E_val(1));
                    subplot(P.Nxplots+length(unicond),P.Nxplots,pln); topoplot(mean(peakdata(:,lat),2), chanlocs,'maplimits','absmax','electrodes','on','plotchans',plotchans);%,'emarker2',{markchans,'o','w',7,1}); 
                    title(num2str(EEGtimes(lat)'))
                end
            end
            %colorbar
        end
        tightfig(f1)
        %save figure 
        if P.save_topo
            print(fullfile(P.spm_path,P.clusdir,[cllabel '_topo']),'-dpng');
        end
    end
        
end
%x=range
end