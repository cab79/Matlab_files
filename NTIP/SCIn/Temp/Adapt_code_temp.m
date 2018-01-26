
% calculate mean of intensity for each type of response in order to
% dynamically adjust the mean over time
if h.Settings.adaptive(atype).getmeanfromresponses
    try
        s.inten_mean = s.out.adaptive(end, 11); % by default, should be same as previous trial, unless modified later
    end
    try
        resp1_ind = find(s.out.adaptive(:, 8)==1); % responses
        resp2_ind = find(s.out.adaptive(:, 8)==2); % responses
        nback = min(h.Settings.adaptive(atype).getmeanfromresponses,size(s.out.adaptive,1));
        ind = length(s.out.adaptive)-nback+1:length(s.out.adaptive);
        resp1_ind = intersect(resp1_ind,ind);
        resp2_ind = intersect(resp2_ind,ind);
        %nback1 = min(length(resp1_ind),h.Settings.adaptive(atype).getmeanfromresponses);
        %nback2 = min(length(resp2_ind),h.Settings.adaptive(atype).getmeanfromresponses);
        % first, calculate the new mean as the mean of recent 1s and 2s
        if ~isempty(resp1_ind) && ~isempty(resp2_ind)
            s.inten_mean(1) = mean(s.out.adaptive(resp1_ind, 9)); % stimulus level
            s.inten_mean(2) = mean(s.out.adaptive(resp2_ind, 9)); % stimulus level
            s.inten_mean = mean(s.inten_mean);
        end
        s.rowofoutput (1, 10) = s.inten_mean; % absolute stimulus intensity
    catch
        s.rowofoutput (1, 10) = NaN;
    end
    
    %if this trial had an incorrect response, evaluate the last few
    %incorrect responses to make mean intenisty adjustments
    if s.rowofoutput (1, 4)==0
        try
            incorr_ind = find(s.out.adaptive(:, 4)==0); % incorrect responses
            resp1_ind = find(s.out.adaptive(:, 8)==1); % responses
            resp2_ind = find(s.out.adaptive(:, 8)==2); % responses
            % second, adjust the mean if ratings are biased in one direction
            nback = min(length(incorr_ind),h.Settings.adaptive(atype).getmeanfromresponses);
            nIncorr1 = length(intersect(resp1_ind,incorr_ind(end-nback+1:end)));
            nIncorr2 = length(intersect(resp2_ind,incorr_ind(end-nback+1:end)));
            % give them a min value of 1 in case there is a strong bias in one
            % direction (meaning the other is never incorrect). This means they
            % have to get more then 1 wrong in a certain direction before
            % adjustments are made.
            if isempty(nIncorr1) || nIncorr1==0
                nIncorr1=1;
            end
            if isempty(nIncorr2) || nIncorr2==0
                nIncorr2=1;
            end
            % adjust in proportion to ratio of incorrect responses
            if nIncorr1<nIncorr2
                s.inten_mean = s.inten_mean-(h.Settings.adaptive(atype).meanadjustmax * (nIncorr2-nIncorr1)/(nIncorr1+nIncorr2));
            elseif nIncorr1>nIncorr2
                s.inten_mean = s.inten_mean+(h.Settings.adaptive(atype).meanadjustmax * (nIncorr1-nIncorr2)/(nIncorr1+nIncorr2));
            end
        end
    end
    try
        s.rowofoutput (1, 11) = s.inten_mean; % absolute stimulus intensity
        disp(['mean adjusted to: ' num2str(s.inten_mean)]);
    catch
        if isfield(s,'inten_mean')
            s = rmfield(s,'inten_mean');
        end
        s.rowofoutput (1, 11) = NaN;
    end
end