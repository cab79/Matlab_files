function factor_img(subdir,imglist,cond_list,scanname,outputtype,norm)
cond = unique(cond_list);
Ncond = length(cond);

switch outputtype
    case 'meancond'
        for c=1:Ncond
            condimg=imglist(cond_list==cond(c));
            fnames={};
            for i = 1:length(condimg)
                fnames{i} = fullfile(subdir, [condimg{i}]); 
                if norm
                    fout = fullfile(subdir,[scanname '_znorm_' num2str(cond(c)) '.nii']);
                    if ~exist(fout,'file')
                        spm_imcalc_ui(fnames,fout,'mean(X)./std(X)',{1,[],[],[]});
                    end
                else
                    fout = fullfile(subdir,[scanname '_' num2str(cond(c)) '.nii']);
                    if ~exist(fout,'file')
                        spm_imcalc_ui(fnames,fout,'mean(X)',{1,[],[],[]});
                    end
                end
            end
        end
    case 'contrast'
        ncondimg = [];
        expr=''
        f=0;
        fnames={};
        for c=Ncond:-1:1
            condimg=imglist(cond_list==cond(c));
            ncondimg(c) = length(condimg);
            for i = 1:ncondimg(c)
                f=f+1;
                % build expression
                if i==1 && c==Ncond
                    expr = '(';
                elseif i==1 && c==1
                    expr = [expr ')/' num2str(ncondimg(c)) '-('];
                else
                    expr = [expr '+'];
                end
                if i==ncondimg(c) && c==1
                    expr = [expr 'i' num2str(f) ')/' num2str(ncondimg(c))];
                else
                    expr = [expr 'i' num2str(f)];
                end
                fnames{f} = fullfile(subdir, [condimg{i}]); 
            end
        end
        fout = fullfile(subdir,[scanname '_contrast.nii']);
        if ~exist(fout,'file')
            spm_imcalc_ui(fnames,fout,expr);
        end
        if norm
            fnorm = fullfile(subdir,'std_allcond.nii');
            if ~exist(fnorm,'file')
                spm_imcalc_ui(fnames,fnorm,'std(X)',{1,[],[],[]});
            end
            fnamesnorm = {fout;fnorm};
            fout = fullfile(subdir,[scanname '_znorm_contrast.nii']);
            if ~exist(fout,'file')
                spm_imcalc_ui(fnamesnorm,fout,'i1./i2');
            end
        end
end