function S=CORE_get_median_params_from_fitted(S,varargin)

temp=load(fullfile(S.path.hgf,'fitted',S.fitted_hgf));
D_fit=temp.D_fit;

if ~isempty(varargin)
    subs = varargin{1};
else
    subs = 1:length(D_fit);
end


% get param fields
sim_fields = fieldnames(D_fit(1).HGF.fit.p_prc); 

% get data and average
for fn = 1:length(sim_fields)
    i=0;
    for d = subs
        i=i+1;
        alldat.(sim_fields{fn})(i,:) = D_fit(d).HGF.fit.p_prc.(sim_fields{fn});
    end
    S.sim.(sim_fields{fn}) = median(alldat.(sim_fields{fn}),1);
end

%test:
%S.sim.like_al0 = S.sim.like_al0 * 2;
%S.sim.PL_om(2) = S.sim.PL_om(2) + 2;
%S.sim.PL_om(3) = S.sim.PL_om(3) - 2;