
    % errors
    out = dir('*.err');
    if ~isempty(out)
        fsize = [out(:).bytes];
        err_ind = find(fsize>0);
        for n = 1:length(err_ind)
            err=load(out(err_ind(n)).name);
            disp(['file size, bytes: ' num2str(out.bytes)])
        end
    end