function fout = factor_img(subdir,imglist,cond_list,scanname,outputtype,norm,trans)

if isempty(cond_list)
    cond_list=1;
end

% condition for summarising
condsum = unique(cond_list(cond_list(:,end)>0,end));
Ncondsum = length(condsum);
if size(cond_list,2)>1
    % condition for contrasting
    condcon = unique(cond_list(cond_list(:,1)>0,1));
    Ncondcon = length(condcon);
end

switch outputtype
    case 'meancond'
        for c=1:Ncondsum
            condimg=imglist(cond_list==condsum(c));
            fnames={};
            for i = 1:length(condimg)
                fname = dir(fullfile(subdir,condimg{i}));
                fnames{i,1} = fullfile(subdir, fname.name); 
            end
            if norm
                fout = fullfile(subdir,[scanname '_znorm_' num2str(condsum(c)) '_' outputtype '.nii']);
                %if ~exist(fout,'file')
                if strcmp(trans,'log')
                    expr = 'log(mean(X)./std(X))'; 
                else
                    expr = 'mean(X)./std(X)'; 
                end
                spm_imcalc_ui(fnames,fout,expr,{1,[],[],[]});
                %end
            else
                fout = fullfile(subdir,[scanname '_' num2str(condsum(c)) '_' outputtype '.nii']);
                %if ~exist(fout,'file')
                if length(fnames)>1
                    if strcmp(trans,'log')
                        expr = 'log(mean(X))'; 
                    else
                        expr = 'mean(X)'; 
                    end
                    spm_imcalc_ui(fnames,fout,expr,{1,[],[],[]});
                elseif strcmp(trans,'log')
                    expr = 'log(i1)'; 
                    spm_imcalc_ui(fnames,fout,expr);
                end
                
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
                % get the images to average for the first (c==2) and second
                % (c==1) part of the contrast
                condimg=imglist(cond_list(:,end)==condsum(c) & cond_list(:,1)==condcon(cs));
                ncondimg(c) = length(condimg);
                for i = 1:ncondimg(c)
                    f=f+1;
                    % build expression
                    if i==1 && c==Ncondsum % if first image of the first part of the contrast
                        if strcmp(trans,'log')
                            expr = 'log('; 
                        else
                            expr = '('; 
                        end
                    elseif i==1 && c==1
                        if strcmp(trans,'log')
                            expr = [expr ')/' num2str(ncondimg(c)) '-log(']; % if the first image of the second part of the contrast
                        else
                            expr = [expr ')/' num2str(ncondimg(c)) '-(']; % if the first image of the second part of the contrast
                        end
                    else
                        expr = [expr '+']; % if second image of any part
                    end
                    if i==ncondimg(c) && c==1
                        expr = [expr 'i' num2str(f) ')/' num2str(ncondimg(c))]; % if second image of second part of contrast
                    else
                        expr = [expr 'i' num2str(f)]; % if second image of first part of contrast
                    end
                    % get file name
                    fname = dir(fullfile(subdir,condimg{i}));
                    fnames{f} = fullfile(subdir, fname.name); 
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