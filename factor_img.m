function factor_img(subdir,imglist,cond_list,scanname)
cond = unique(cond_list);
Ncond = length(cond);
for c=1:Ncond
    condimg=imglist(cond_list==cond(c));
    fnames={};
    for i = 1:length(condimg)
        fnames{i} = fullfile(subdir, [condimg{i}]);
    end
    fout = fullfile(subdir,[scanname '_' num2str(cond(c)) '.nii']);
    spm_imcalc_ui(fnames,fout,'mean(X)',{1,[],[],[]});
end