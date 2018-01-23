function fout = factor_img(subdir,imglist,cond_list,scanname,outputtype,norm)
% condition for summarising
condsum = unique(cond_list(cond_list(:,end)>0));
Ncondsum = length(condsum);
if size(cond_list,2)>1
    % condition for contrasting
    condcon = unique(cond_list(cond_list(:,1)>0));
    Ncondcon = length(condcon);
end

switch outputtype
    case 'meancond'
        for c=1:Ncondsum
            condimg=imglist(cond_list==condsum(c));
            fnames={};
            for i = 1:length(condimg)
                fnames{i,1} = fullfile(subdir, [condimg{i}]); 
            end
            if norm
                fout = fullfile(subdir,[scanname '_znorm_' num2str(condsum(c)) '_' outputtype '.nii']);
                if ~exist(fout,'file')
                    spm_imcalc_ui(fnames,fout,'mean(X)./std(X)',{1,[],[],[]});
                end
            else
                fout = fullfile(subdir,[scanname '_' num2str(condsum(c)) '_' outputtype '.nii']);
                %if ~exist(fout,'file')
                    spm_imcalc_ui(fnames,fout,'mean(X)',{1,[],[],[]});
                %end
            end
        end
    case 'meanallcond'
            condimg=imglist;
            fnames={};
            for i = 1:length(condimg)
                fnames{i,1} = fullfile(subdir, [condimg{i}]); 
            end
            if norm
                fout = fullfile(subdir,[scanname '_znorm_' num2str(condsum(c)) '_' outputtype '.nii']);
                if ~exist(fout,'file')
                    spm_imcalc_ui(fnames,fout,'mean(X)./std(X)',{1,[],[],[]});
                end
            else
                fout = fullfile(subdir,[scanname '_' outputtype '.nii']);
                %if ~exist(fout,'file')
                    spm_imcalc_ui(fnames,fout,'mean(X)',{1,[],[],[]});
                %end
            end
    case 'contrast'
        ncondimg = [];
        expr=''
        f=0;
        fnames={};
        for cs = 1:Ncondcon
            for c=Ncondsum:-1:1
                condimg=imglist(cond_list(:,end)==condsum(c) & cond_list(:,1)==condcon(cs));
                ncondimg(c) = length(condimg);
                for i = 1:ncondimg(c)
                    f=f+1;
                    % build expression
                    if i==1 && c==Ncondsum
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
                    fnames{f} = fullfile(subdir, condimg{i}); 
                end
            end
            fout = fullfile(subdir,[scanname num2str(cs) '_' outputtype '.nii']);
            %if ~exist(fout,'file')
                spm_imcalc_ui(fnames,fout,expr);
            %end
            if norm
                fnorm = fullfile(subdir,'std_allcond.nii');
                if ~exist(fnorm,'file')
                    spm_imcalc_ui(fnames,fnorm,'std(X)',{1,[],[],[]});
                end
                fnamesnorm = {fout;fnorm};
                fout = fullfile(subdir,[scanname '_znorm_' outputtype '.nii']);
                if ~exist(fout,'file')
                    spm_imcalc_ui(fnamesnorm,fout,'i1./i2');
                end
            end
        end
end