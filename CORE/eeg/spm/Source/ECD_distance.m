% digit_mni_prior = [
%     44, -22.8, 59.9;
%     40.6, -28.2, 62.1;
%     37.7, -29.2, 64.8;
%     35.4, -30, 66.3;
%     ];


restoredefaultpath
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')
dbstop if error
outpath = 'C:\Data\CORE\eeg\ana\spm\SPMdata\ECD_outputs';
cd(outpath)
set(0,'DefaultFigureVisible','on')


%file = 'ECD_out_20181108T192602.mat';
%file = 'ECD_out_20181110T073609.mat'; % 39.4  -27.6   63.3 (mean), 
%file = 'ECD_out_20181110T065516.mat'; % 24 -34 64  (mean), 
%file = 'ECD_out_20181110T233644.mat'; % peak 1
%file = 'ECD_out_20181110T234844.mat'; % peak 2
% 3 dipoles
%file = 'ECD_out_20181113T050045.mat'; % peak 1
%file = 'ECD_out_20181113T013209.mat'; % peak 2
file = 'ECD_out_20181113T020012.mat'; % peak 3

load(file);

hand = {[1:4],[5:8]};
dipole = 1; % index of dipole to analyse

% re-arrange data for stats and plots
dist=struct;
stats=struct;
loc=struct;
for g = 1:2 % groups
    nsub = size(mniloc{g},1);
    locmom{g} = cat(3,mniloc{g},jmni{g});
    
    for h = 1:2 % hand
        for sub = 1:nsub
            X = squeeze(mniloc{g}(sub,hand{h},:,dipole));
            for d = 1:size(X,1) % digit pair
                loc(h,d).grp{g}(sub,:) = X(d,:);
            end

            loc_dist = pdist(X);
            for d = 1:length(loc_dist) % digit pair
                dist(h,d).loc{g}(sub,1) = loc_dist(d);
            end

            X = squeeze(jmni{g}(sub,hand{h},:));
            mom_dist = pdist(X);
            for d = 1:length(mom_dist) % digit pair
                dist(h,d).mom{g}(sub,1) = mom_dist(d);
            end
        end
    end
end
grps = ([ones(size(loc(1,1).grp{1},1),1); 2*ones(size(loc(1,1).grp{2},1),1)]);

% plot
close all
plot_var=0; % set to 0 to plot digits in different colours
for g = 1:2 % groups
    for h = 1:2 % hand
        figure
        
        % prior
        mloc = S.prior(dipole).prior_loc;
        vloc = S.prior(dipole).prior_var/2;
        if plot_var
            [x, y, z] = ellipsoid(mloc(1),mloc(2),mloc(3),vloc(1),vloc(2),vloc(3),30);
            surf(x, y, z,'FaceAlpha',0.2,'EdgeColor','none','FaceColor','k')
            hold on
            scatter3(mloc(1),mloc(2),mloc(3),100,'k','filled');
        end
        %axis equal
        
        % estimated dipoles
        x=[];y=[];z=[];
        for d=1:4
            x = [x;loc(h,d).grp{g}(:,1)];
            y = [y;loc(h,d).grp{g}(:,2)];
            z = [z;loc(h,d).grp{g}(:,3)];
        end
        % shape and colour
        Sh = repmat([50,50,50,50],numel(x)/d,1);
        C = repmat([1,2,3,4],numel(x)/d,1);
        s = Sh(:);
        c = C(:);
        % plot
        scatter3(x,y,z,s,c);
        %hnd.MarkerFaceColor = [0 0.5 0.5];
        title(['dipole ' num2str(dipole) ', group '  num2str(g) ', hand ' num2str(h)])
        
     end
end
drawnow
%pause(1)

% save on T1 image
node_radius = 3;
V=spm_vol('C:\Data\Matlab\spm12\canonical\single_subj_T1.nii');
img=spm_read_vols(V);
cimg = 1:numel(img);
[ix iy iz] = ind2sub(size(img),cimg');
ind3 = [ix iy iz];
mni = cor2mni(ind3, V.mat); % convert to MNI
for g = 1:2 % groups
    for h = 1:2 % hand
        for d=1:4
            save_img = zeros(size(img));
            subj_nodes = loc(h,d).grp{g};
            dNodes2Centers = pdist2(mni, subj_nodes);
            [ro, co] = find(dNodes2Centers < node_radius);
            nodeIdsPerCenter = accumarray(co, ro, [], @(x){x});
            for sub = 1:length(nodeIdsPerCenter)
                save_img(nodeIdsPerCenter{sub}) = d;
            end
            Vs=V;
            fname = strrep(file,'.mat',['_dip' num2str(dipole) '_grp' num2str(g) '_hand' num2str(h) '_digit' num2str(d) '.nii']);
            Vs.fname=fullfile(outpath,fname);
            Vs.pinfo=[1;0;0];
            spm_write_vol(Vs,save_img);
        end
    end
end

% stats: Digit
digit_pairs = nchoosek(1:4,2);
for h=1:2
    for d=1:size(digit_pairs,1)
        % both groups
        dat = [loc(h,digit_pairs(d,1)).grp{1};loc(h,digit_pairs(d,1)).grp{2};loc(h,digit_pairs(d,2)).grp{1};loc(h,digit_pairs(d,2)).grp{2}];
        conds = [ones(size(grps,1),1); 2*ones(size(grps,1),1)];
        out = lda_class(dat,conds,length(conds),5,100);
        stats.loc.cond.lda(h,d) = out.ldaCVErr;
        stats.loc.cond.lda_pval(h,d) = out.pval;
        % CRPS
        dat = [loc(h,digit_pairs(d,1)).grp{1};loc(h,digit_pairs(d,2)).grp{1}];
        conds = [ones(size(loc(1,1).grp{1},1),1); 2*ones(size(loc(1,1).grp{1},1),1)];
        out = lda_class(dat,conds,length(conds),5,100);
        stats.loc.cond_grp1.lda(h,d) = out.ldaCVErr;
        stats.loc.cond_grp1.lda_pval(h,d) = out.pval;
        % HC
        dat = [loc(h,digit_pairs(d,1)).grp{2};loc(h,digit_pairs(d,2)).grp{2}];
        conds = [ones(size(loc(1,1).grp{2},1),1); 2*ones(size(loc(1,1).grp{2},1),1)];
        out = lda_class(dat,conds,length(conds),5,100);
        stats.loc.cond_grp2.lda(h,d) = out.ldaCVErr;
        stats.loc.cond_grp2.lda_pval(h,d) = out.pval;
    end
end

% stats: Grp
%grps = {'CRPS','HC'};
for h=1:2
    for d=1:4
        stats.loc.ranksum(h,d) = ranksum(dist(h,d).loc{1},dist(h,d).loc{2});
        stats.mom.ranksum(h,d) = ranksum(dist(h,d).mom{1},dist(h,d).mom{2});
        out = lda_class([loc(h,d).grp{1};loc(h,d).grp{2}],grps,length(grps),5,100);
        stats.loc.grp.lda(h,d) = out.ldaCVErr;
        stats.loc.grp.lda_pval(h,d) = out.pval;
    end
end
