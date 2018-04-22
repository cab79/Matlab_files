function h = trial_set_param(h)

switch h.Settings.stim(h.trialstimnum).control
    
    case {'LJTick-DAQ','labjack'}

        if ~isfield(h.Settings.stim,'inten_diff')
            h.Settings.stim(h.trialstimnum).inten_diff = 0;
        end
        h.dur = h.Settings.stim(h.trialstimnum).dur; 
       
        % set initial values from GUI or settings if not already
        % defined
        set_inten_mean = 0;
        if ~isfield(h,'inten_mean')
            set_inten_mean = 1;
        elseif isempty(h.inten_mean)
            set_inten_mean = 1;
        end
        if set_inten_mean
            if ~isfield(h,'inten_mean_gui')
                h.inten_mean=0;
            else
                h.inten_mean = str2double(h.inten_mean_gui);
            end
            if ~any(h.inten_mean) && ~isempty(h.Settings.stim(h.trialstimnum).inten)
                h.inten_mean = h.Settings.stim(h.trialstimnum).inten;
            end
            h.inten_mean_start = h.inten_mean;
        end
        set_inten_diff = 0;
        if ~isfield(h,'inten_diff')
            set_inten_diff = 1;
        elseif isempty(h.inten_diff)
            set_inten_diff = 1;
        end
        if set_inten_diff
            if ~isfield(h,'inten_diff_gui')
                h.inten_diff=0;
            else
                h.inten_diff = str2double(h.inten_diff_gui);
            end
            if ~any(h.inten_diff) && ~isempty(h.Settings.stim(h.trialstimnum).inten_diff)
                h.inten_diff = h.Settings.stim(h.trialstimnum).inten_diff;
            end
            h.inten_diff_start = h.inten_diff;
        end

        % modify according to procedure
        if h.seqtype.thresh
            if ~isfield(h,'s')
                h.inten = h.inten_mean+h.Settings.threshold.startinglevel;
            else
                h.inten = h.inten_mean+h.s.StimulusLevel;
            end
        % adaptive trial
        elseif h.seqtype.adapt
            detect = find(strcmp(h.atypes,'detect'));
            discrim = find(strcmp(h.atypes,'discrim'));
            % if this trial is adaptive
            if ~isnan(h.Seq.adapttype(h.i))
                % update mean from adaptive procedure.
                % do this even if it's a discrim trial
                if ~isempty(detect) && isfield(h,'s') 
                    if isfield(h.s.a(detect),'StimulusLevel') % if a level has been determined
                        h.inten_mean = h.s.a(detect).StimulusLevel;
                        if isempty(h.inten_mean)
                            h.inten_mean = h.inten_mean_start;
                        end
                    end
                % or set the adaptive starting level otherwise
                elseif ~isfield(h,'s')
                    h.Settings.adaptive(detect).startinglevel = h.inten_mean;
                end

                % update diff from adaptive
                if ~isempty(discrim)
                    if isfield(h,'s') && length(h.s.a)>=discrim % if a level has been determined
                        h.inten_diff = h.s.a(discrim).StimulusLevel;
                    end
                    if h.inten_diff == 0 || isempty(h.inten_diff) % if not set in GUI or settings
                        h.inten_diff = h.Settings.adaptive(discrim).startinglevel;
                    end
                end

                % only do this if it's a discrim trial, not a detect trial
                if h.Seq.adapttype(h.i) == discrim
                    % calculate intensity
                    if h.Seq.signal(h.trialstimnum,h.i)==h.trialstimnum
                        h.inten = h.inten_mean + h.Settings.adaptive(discrim).stepdir * h.inten_diff / 2;
                    else
                        h.inten = h.inten_mean - h.Settings.adaptive(discrim).stepdir * h.inten_diff / 2;
                    end
                else
                    h.inten = h.inten_mean;
                end

            % if adaptive is part of the sequence, but not this trial
            elseif isnan(h.Seq.adapttype(h.i))
                detect_thresh =  find(h.out.adaptive(:,10)==detect);
                discrim_thresh =  find(h.out.adaptive(:,10)==discrim);
                if h.Seq.signal(h.trialstimnum,h.i)==1
                    h.inten = h.out.adaptive(detect_thresh(end),7) - h.out.adaptive(discrim_thresh(end),7) / 2;
                else
                    h.inten = h.out.adaptive(detect_thresh(end),7) + h.out.adaptive(discrim_thresh(end),7) / 2;
                end
            end
        % Otherwise, use sequence to determine intensity
        else
            if strcmp(h.Settings.oddballmethod,'intensity')
                if iscell(h.Settings.oddballvalue)
                    if size(h.Settings.oddballvalue,1)==1
                        h.inten = h.Settings.oddballvalue{h.Seq.signal(h.trialstimnum,h.tr)};
                    else
                        h.inten = h.Settings.oddballvalue{h.Seq.signal(h.trialstimnum,h.tr),:};
                    end
                else
                    h.inten = h.Settings.oddballvalue(h.Seq.signal(h.trialstimnum,h.tr),:);
                end
            elseif strcmp(h.Settings.oddballmethod,'index')
                % calculate intensity
                if h.Seq.signal(h.trialstimnum,h.i)==1
                    h.inten = h.inten_mean - h.inten_diff / 2;
                else
                    h.inten = h.inten_mean + h.inten_diff / 2;
                end
            end
        end

        % set max intensity
        if ~isfield(h.Settings.stim,'maxinten')
            h.Settings.stim(h.trialstimnum).maxinten = inf;
        end
        h.inten = min(h.inten,h.Settings.stim(h.trialstimnum).maxinten); 
        h.inten_out = h.inten;

    case {'PsychPortAudio','audioplayer'}
    
        % create temporary variables for intensity, pitch and duration
        h.inten = h.Settings.stim(h.trialstimnum).inten;
        h.freq = h.Settings.stim(h.trialstimnum).f0;
        h.dur = h.Settings.stim(h.trialstimnum).dur; 
        h.varlevel = 0;
        % if calculating all trials, requires totdur:
        if strcmp(h.Settings.design,'continuous') && ~isfield(h,'totdur') % Use calculated duration by default.
            if isfield(h.Settings,'totdur')
                h.totdur = h.Settings.totdur; 
            end
        end
        if ~isfield(h.Settings.stim(h.trialstimnum),'inten_type')
            h.Settings.stim(h.trialstimnum).inten_type = 'abs';
        end
        if ~isfield(h,'inten_diff')
            if ~isfield(h,'aud_diff_gui')
                h.inten_diff=0;
            else
                h.inten_diff = str2double(h.aud_diff_gui);
            end
            if ~any(h.inten_diff) && ~isempty(h.Settings.stim(h.trialstimnum).inten_diff)
                h.inten_diff = h.Settings.stim(h.trialstimnum).inten_diff;
            end
            h.inten_diff_start = h.inten_diff;
        end
        
        % thresholding
        if h.seqtype.thresh && ~h.seqtype.oddball
            if ~isfield(h,'s')
                h.varlevel = h.Settings.threshold.startinglevel;
            else
                h.varlevel = h.s.StimulusLevel;
            end
        end

        % condition method
        if isfield(h.Settings,'conditionmethod')
            if ~isempty(h.Settings.conditionmethod)
                if iscell(h.Settings.conditionmethod)
                    for i = 1:length(h.Settings.conditionmethod)
                        conditionmethod = h.Settings.conditionmethod{i};
                        if strcmp(conditionmethod,'pitch') || strcmp(conditionmethod,'freq')
                            if iscell(h.Settings.conditionvalue)
                                h.freq = h.Settings.conditionvalue{h.Seq.signal(h.trialstimnum,h.tr),i};
                            else
                                h.freq = h.Settings.conditionvalue(i,h.Seq.signal(h.trialstimnum,h.tr));
                            end
                        end
                        if strcmp(conditionmethod,'intensity')
                            if iscell(h.Settings.conditionvalue)
                                h.inten = h.Settings.conditionvalue{h.Seq.signal(h.trialstimnum,h.tr),i};
                            else
                                h.inten = h.Settings.conditionvalue(i,h.Seq.signal(h.trialstimnum,h.tr));
                            end
                        end
                        if strcmp(conditionmethod,'phase')
                            if iscell(h.Settings.conditionvalue)
                                h.alignphase = h.Settings.conditionvalue{h.Seq.signal(h.trialstimnum,h.tr),i};
                            else
                                h.alignphase = h.Settings.conditionvalue(i,h.Seq.signal(h.trialstimnum,h.tr));
                            end
                        end
                    end
                else
                    error('h.Settings.conditionmethod must be a cell');
                end
            end
        end

        % oddball method
        if h.seqtype.oddball
            if ~h.seqtype.adapt && ~h.seqtype.thresh
                if iscell(h.Settings.oddballvalue)
                    if size(h.Settings.oddballvalue,1)==1
                        oddval = h.Settings.oddballvalue{h.Seq.signal(h.trialstimnum,h.tr)};
                    else
                        oddval = h.Settings.oddballvalue{h.Seq.signal(h.trialstimnum,h.tr),:};
                    end
                else
                    oddval = h.Settings.oddballvalue(h.Seq.signal(h.trialstimnum,h.tr),:);
                end
                if strcmp(h.Settings.oddballmethod,'channel')
                    h.chan=oddval;
                elseif strcmp(h.Settings.oddballmethod,'intensity') && strcmp(h.Settings.stim(h.trialstimnum).inten_type,'abs') && h.inten_diff==0
                    h.inten = oddval;
                elseif strcmp(h.Settings.oddballmethod,'index') && strcmp(h.Settings.stim(h.trialstimnum).inten_type,'dB') && h.inten_diff~=0 && strcmp(h.Settings.stim(h.trialstimnum).patternmethod,'intensity') && isempty(h.Settings.stim(h.trialstimnum).patternvalue)
                    if isfield(h.Settings.assoc,'intenstim')
                        h.varlevel = h.Settings.assoc.intenstim(h.Seq.signal(h.trialstimnum,h.i)) * h.inten_diff/2;
                    else
                        if h.Seq.signal(h.trialstimnum,h.i)==1
                            h.varlevel = h.inten_diff/2;
                        else
                            h.varlevel = -h.inten_diff/2;
                        end
                    end
                elseif strcmp(h.Settings.oddballmethod,'pitch') || strcmp(h.Settings.oddballmethod,'freq')
                    h.freq = oddval;
                elseif strcmp(h.Settings.oddballmethod,'duration')
                    h.dur = oddval;
                end
            elseif (h.seqtype.adapt && h.Settings.adaptive_general.stim==h.trialstimnum) || h.seqtype.thresh
                if isfield(h,'s')
                    h.varlevel = h.s.a(h.Seq.adapttype(h.i)).StimulusLevel;
                else
                    if h.seqtype.adapt
                        h.varlevel = h.Settings.adaptive.startinglevel;
                    else
                        h.varlevel = h.Settings.threshold.startinglevel;
                    end  
                end
                if strcmp(h.Settings.oddballmethod,'pitch') || strcmp(h.Settings.oddballmethod,'freq')
                    h.freq = [h.Settings.oddballvalue(1), (h.Settings.oddballvalue(1)+h.varlevel)]; % create new pitch pair
                    h.freq = h.freq(h.Seq.signal(h.trialstimnum,h.tr));
                elseif strcmp(h.Settings.oddballmethod,'intensity') && strcmp(h.Settings.stim(h.trialstimnum).inten_type,'abs')
                    h.inten = [h.Settings.oddballvalue(1), (h.Settings.oddballvalue(1)+h.varlevel)]; % create new pitch pair
                    h.inten = h.inten(h.Seq.signal(h.trialstimnum,h.tr));
                elseif strcmp(h.Settings.oddballmethod,'duration') && (strcmp(h.Settings.stim(h.trialstimnum).patternmethod,'pitch') || strcmp(h.Settings.stim(h.trialstimnum).patternmethod,'freq'))
                    if iscell(h.Settings.oddballvalue)
                        h.dur = h.Settings.oddballvalue{h.Seq.signal(h.trialstimnum,h.tr),:};
                    else
                        h.dur = h.Settings.oddballvalue(h.Seq.signal(h.trialstimnum,h.tr),:);
                    end
                    h.freq = [h.Settings.stim(h.trialstimnum).patternvalue(1), (h.Settings.stim(h.trialstimnum).patternvalue(1)+h.varlevel)]; % create new pitch pair
                end
            end
        end

        %apply pitch pattern?
        h.trialtype.freqpattern=0;
        if isfield(h.Settings.stim(h.trialstimnum),'patternmethod')
            if strcmp(h.Settings.stim(h.trialstimnum).patternmethod,'pitch') || strcmp(h.Settings.stim(h.trialstimnum).patternmethod,'freq') % pitch changes
                h.trialtype.freqpattern=1;
                if ~((h.seqtype.adapt || h.seqtype.thresh) && (strcmp(h.Settings.oddballmethod,'pitch') || strcmp(h.Settings.oddballmethod,'freq'))) && ~(~isempty(strcmp(h.Settings.conditionmethod,'pitch')) || ~isempty(strcmp(h.Settings.conditionmethod,'freq'))) % pitch already defined above in this case
                    if isnumeric(h.Settings.stim(h.trialstimnum).patternvalue)
                        h.freq = h.Settings.stim(h.trialstimnum).patternvalue;
                    elseif iscell(h.Settings.stim(h.trialstimnum).patternvalue) && strcmp(h.Settings.oddballmethod,'index')
                        h.freq = h.Settings.stim(h.trialstimnum).patternvalue{h.Seq.signal(h.trialstimnum,h.i)};
                    end
                    % OLD CODE: NOT NEEDED?
                    %elseif iscell(h.Settings.stim(h.trialstimnum).patternvalue)
                    %    nDur = length(h.dur);
                    %    nPit = cellfun(@length,h.Settings.stim(h.trialstimnum).patternvalue);
                    %    h.freq = h.Settings.stim(h.trialstimnum).patternvalue{nPit==nDur};
                    %end
                end
            end
        end
        % apply response probe?
        if isfield(h.Settings,'RPmethod')
            if strcmp(h.Settings.RPmethod,'pitch') || strcmp(h.Settings.RPmethod,'freq')
                if h.Seq.RP(h.tr)==1
                    h.trialtype.freqpattern=1;
                    h.freq = h.Settings.RPvalue;
                    h.dur = h.Settings.RPdur;
                end
            end
        end
        %apply intensity pattern?
        h.trialtype.intenpattern=0;
        if isfield(h.Settings.stim(h.trialstimnum),'patternmethod')
            if strcmp(h.Settings.stim(h.trialstimnum).patternmethod,'intensity') && ~isempty(h.Settings.stim(h.trialstimnum).patternvalue) % intensity changes
                h.trialtype.intenpattern=1;
                if ~any(strcmp(h.Settings.conditionmethod,'intensity')) % then already defined
                    h.inten = h.Settings.stim(h.trialstimnum).patternvalue;
                end
            end
        end
        % apply response probe?
        h.resp_probe=0;
        if isfield(h.Settings,'RPmethod')
            if strcmp(h.Settings.RPmethod,'intensity')
                if h.Seq.RP(h.tr)==1
                    h.resp_probe=1;
                    h.inten = h.Settings.RPvalue;
                    h.dur = h.Settings.RPdur;
                end
            end
        end
        
end