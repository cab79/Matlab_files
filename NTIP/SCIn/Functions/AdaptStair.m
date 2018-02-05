function h = AdaptStair(h,varargin)

% Each type will calculate their respective thresholds/levels independently
% by using information across blocks of the same type, and will use values from the blocks of the other type.

adaptive=0;
opt='';

% find out whether there is more than one type of adaptive in this sequence
if nargin>1
    atype = varargin{1};
else
    atype = 1;
end

% create s.a?
create_a = 1;
create_s = 1;
if isfield(h,'s')
    create_s = 0;
    if isfield(h.s,'a')
        if length(h.s.a)==atype
            create_a = 0;
        end
    end
end
if create_s || create_a
    s = AdaptStairParameters(h,atype); % see function at the end
else
    s=h.s;
end

%%%%%%%%%%%%%%%%% TEMP - need to sort these out
s.block=1;
%%%%%%%%%%%%%%%%%%

if ~isfield(s,'lasttrialrun')
    s.lasttrialrun = 0;
end

s.mintrialcount = h.i;
 % identify BUTTON PRESSES or OMISSIONS over previous trials of
% the same type
if isfield(h.Settings.adaptive(atype),'oddonly')
    if h.Settings.adaptive(atype).oddonly

        % only continue if next trial (that can be modified! i.e. h.i+h.Settings.ntrialsahead-1) is an oddball
        try
            if h.i+h.Settings.ntrialsahead-1 < length(h.Seq.signal)
                all_odd_ind = horzcat(h.Seq.odd_ind{:});
                %if any(ismember(all_odd_ind,h.Seq.condnum(h.i+h.Settings.ntrialsahead-1+1)))
                %    disp(['next (modifiable) trial is an oddball: trial ' num2str(h.i+h.Settings.ntrialsahead-1+1) ', condnum ' num2str(h.Seq.condnum(h.i+h.Settings.ntrialsahead-1+1))])
                %end
                if h.Settings.adaptive(atype).oddonly && ~any(ismember(all_odd_ind,h.Seq.condnum(h.i+h.Settings.ntrialsahead-1+1))) % h.odd_ind is generated in CreateSequence
                    return
                end
            else
                return
            end
        catch
            error('define h.Seq.odd_ind and h.Seq.condnum in Sequence creation function');
        end

        % any responses over previous trials of this stimtype?
        notthis = find(h.Seq.signal~=h.Seq.signal(h.i)); % trials that are not this stimtype
        prevnotthis = notthis(notthis<h.i);% previous trials that were not this stimtype
        lastnotthis = max(prevnotthis); % last trial
        trls = lastnotthis+1:max(lastnotthis+h.Settings.adaptive(atype).resptrials,h.i); % count to h.i or the max number of allowed trials
        trls_pressed = h.out.presstrial(ismember(h.out.presstrial,trls));
        % if there has been a button press that has not be
        % adapted to yet:
        if ~isempty(trls_pressed)
            h.pressedsinceoddball = 1;
        else
            h.omittedsinceoddball = 1;
        end

        % number of oddballs so far
        s.mintrialcount = sum(ismember(h.Seq.condnum(1:h.i),all_odd_ind));
        
        % run Adaptive if there have been BUTTON PRESSES on
        % this trial or over previous trials that not yet been adapted to
        if h.pressedsinceoddball && s.lasttrialrun~=trls_pressed(end)
            opt = 'responded';
            adaptive=1;
        end
    end
end

if (h.pressed || h.pressedlasttrial || h.pressedthistrial) && s.lasttrialrun~=h.i
    opt = 'responded';
    adaptive=1;
end

% if stimulation requires adaptive tuning to OMISSION of responses
if h.Settings.adaptive(atype).omit && s.mintrialcount>h.Settings.adaptive(atype).startomit
    if h.Settings.adaptive(atype).oddonly
        if h.omittedsinceoddball
            opt = 'omitted';
            adaptive=1;
        end
    elseif h.omittedlasttrial && s.lasttrialrun~=h.i
        opt = 'omitted';
        adaptive=1;
    end
end

if ~adaptive
    return
end

% adapt to omitted response?
%omit=0;
%if ~isempty(varargin)
%   if strcmp(varargin{1},'omitted');
%       omit=1;
%   end
%end

%% RUN
disp(['Running Adaptive: ' h.Settings.adaptive(atype).type])

% only run once per trial
%if strcmp(opt,'responded')
    if ~(s.lasttrialrun == h.i)
        s.lasttrialrun = h.i;
    else
        return
    end
%else
%    s.lasttrialrun = h.i;
%end

if size(s.a(atype).expplan,1)==s.a(atype).count_of_n_of_reversals
    s.a(atype).expplan(end+1,:) = s.a(atype).expplan(end,:);
end
    
% evaluate the subject's response
if strcmp(opt,'omitted')
    s.SubjectAccuracy(s.trial)= 0;
elseif strcmp(opt,'responded')
    if h.pressedsinceoddball 
        %presstrial=trls_pressed(end);
        pressbutton = h.out.pressbutton(h.out.presstrial==trls_pressed(end));
        correctsignal = h.Seq.signal(trls_pressed(end));
    elseif h.pressedlasttrial
        %presstrial=h.i;
        pressbutton = h.out.lastpressbutton{h.i};
        correctsignal = h.Seq.signal(h.i);
    elseif h.pressedthistrial
        %presstrial=h.i;
        pressbutton = h.out.pressbutton(h.out.presstrial==h.i);
        correctsignal = h.Seq.signal(h.i);
    elseif h.pressed
        %presstrial=h.i;
        pressbutton = h.out.pressbutton(h.out.presstrial==h.i);
        correctsignal = h.Seq.signal(h.i);
    end
    if ~iscell(pressbutton)
        pressbutton = {pressbutton};
    end
    if isempty(pressbutton)
        return
    end
    resi = find(strcmp(pressbutton(end),h.Settings.buttonopt)); % which button was pressed?
    resfun = h.Settings.adaptive(atype).signalval(resi); %what is meaning of this response?
    if isempty(resi)
        return
    end
    if resfun == correctsignal && strcmp(h.Settings.adaptive(atype).type,'discrim')
    %s.SubjectAccuracy(s.trial)= EvaluateAnswer(CorrectAnswer,s.feedback,Question);   %evaluing the subject answer (right or wrong)
        s.SubjectAccuracy(s.trial)= 1;
    else
        s.SubjectAccuracy(s.trial)= 0;
    end
end

% UPDATE THE ROWOFOUTPUT
s.rowofoutput (1, 1) = s.block;
s.rowofoutput (1, 2) = s.trial;
s.rowofoutput (1, 3) = s.a(atype).StimulusLevel;
s.rowofoutput (1, 4) = s.SubjectAccuracy(s.trial);
% update the count for the up-down motion
go_down=0;go_up=0;
if strcmp(h.Settings.adaptive(atype).type,'discrim') && s.SubjectAccuracy(s.trial) == 1
    go_down=1;
elseif strcmp(h.Settings.adaptive(atype).type,'detect') && resfun == 2
    go_down=1;
else
    go_up=1;
end
if go_down
    s.a(atype).n_down = s.a(atype).n_down + 1;
    if s.a(atype).n_down == s.p(atype).down
        s.a(atype).n_down = 0;
        s.a(atype).pos = 1;
        s.a(atype).trend = 1;
        % update the count of the number of s.reversals and
        % corresponding stepsize
        if s.a(atype).pos ==1 && s.a(atype).neg == -1
            s.a(atype).count_of_n_of_reversals = s.a(atype).count_of_n_of_reversals + 1;
            % calculate the threshold
            s.a(atype).blockthresholds(s.a(atype).n_threshold)=s.a(atype).StimulusLevel;
            s.a(atype).n_threshold = s.a(atype).n_threshold + 1;
            s.a(atype).actualstep = s.a(atype).expplan(s.a(atype).count_of_n_of_reversals, 2);
            s.a(atype).pos = s.a(atype).trend;
            s.a(atype).neg = s.a(atype).trend;
        end
        if s.p(atype).isstep == 1
            s.a(atype).StimulusLevel = s.a(atype).StimulusLevel - s.a(atype).actualstep;
        else
            s.a(atype).StimulusLevel = s.a(atype).StimulusLevel / s.a(atype).actualstep;
        end
    end
elseif go_up
    %error(['stopped at mintrialcount: ' num2str(s.mintrialcount)])
    s.a(atype).neg = -1;
    s.a(atype).trend = -1;
    s.a(atype).n_down = 0;
    % update the count of the number of s.reversals and
    % corresponding stepsize
    if s.a(atype).pos ==1 && s.a(atype).neg == -1
        s.a(atype).count_of_n_of_reversals = s.a(atype).count_of_n_of_reversals + 1;
        % calculate the threshold
        s.a(atype).blockthresholds(s.a(atype).n_threshold)=s.a(atype).StimulusLevel;
        s.a(atype).n_threshold = s.a(atype).n_threshold + 1;
        s.a(atype).actualstep = s.a(atype).expplan(s.a(atype).count_of_n_of_reversals, 2);
        s.a(atype).pos = s.a(atype).trend;
        s.a(atype).neg = s.a(atype).trend;
    end
    if s.p(atype).isstep == 1
        s.a(atype).StimulusLevel = s.a(atype).StimulusLevel + s.a(atype).actualstep;
    else
        s.a(atype).StimulusLevel = s.a(atype).StimulusLevel * s.a(atype).actualstep;
    end
    if isfield(h.Settings.adaptive(atype),'levelmax')
        s.a(atype).StimulusLevel = min(s.a(atype).StimulusLevel,h.Settings.adaptive(atype).levelmax);
    end
end

% UPDATE THE ROWOFOUTPUT
s.rowofoutput (1, 5) = s.a(atype).count_of_n_of_reversals;
s.rowofoutput (1, 6) = s.a(atype).actualstep;
% update the number of trials
s.trial = s.trial + 1;

%disp(['length_blockthresh = ' num2str(length(s.blockthresholds))]);
disp(['nreversals = ' num2str(s.a(atype).count_of_n_of_reversals)]);
disp(['next level = ' num2str(s.a(atype).StimulusLevel)]);

% threshold for the block
if length(s.a(atype).blockthresholds)>=s.p(atype).reversalForthresh
    switch s.p(atype).thresholdtype
        case 'Arithmetic'
            s.a(atype).expthresholds(s.block)=mean(s.a(atype).blockthresholds(end-(s.p(atype).reversalForthresh-1):end));
        case 'Geometric'
            s.a(atype).expthresholds(s.block)=prod(s.a(atype).blockthresholds(end-(s.p(atype).reversalForthresh-1):end))^(1/length(s.a(atype).blockthresholds(end-(s.p(atype).reversalForthresh-1):end)));
        case 'Median'
            s.a(atype).expthresholds(s.block)=median(s.a(atype).blockthresholds(end-(s.p(atype).reversalForthresh-1):end));
        otherwise
            disp('Unknown calculation type.')
    end
    fprintf('Threshold equal to %1.3f\n', s.a(atype).expthresholds(s.block));
else
    s.a(atype).expthresholds(s.block)=NaN;
end
s.rowofoutput (1, 7) = s.a(atype).expthresholds(s.block);
s.rowofoutput (1, 8) = resfun; % the actual response meaning
s.rowofoutput (1, 9) = h.inten; % absolute stimulus intensity
s.rowofoutput (1, 10) = atype; % absolute stimulus intensity


% UPDATE THE GLOBAL OUTPUT VARIABLE
s.out.adaptive = [s.out.adaptive; s.rowofoutput];
h.s =s;
h.out.adaptive = s.out.adaptive;
%fprintf('Press return to continue\n');
%pause
%fprintf ('\nBLOCK ENDED\n');
%pause(2)

%% plot
ind = find(h.out.adaptive(:,10)==atype & ~isnan(h.out.adaptive(:,7)));
if ~isempty(ind)
    select_ind = ind(end-(min(h.Settings.adaptive(atype).av_thresh,length(ind)))+1,end);
    thresh = h.out.adaptive(ind,7);
    av_thresh = mean(thresh);

    % does fig handle exist?
    fig = isfield(h,'f');
    if fig
        fig = fig && length(h.f)>=atype;
    end
    if fig
        eval(['fig = ishandle(h.f(' num2str(atype) '));']);
    end
    
    if ~fig
        set(groot, 'DefaultFigureVisible', 'off');
        eval(['h.f(' num2str(atype) ')=figure;']);
    else
        eval(['set(groot, ''CurrentFigure'', h.f(' num2str(atype) '));']);
        %eval(['figure(h.f' num2str(atype) ');']);
    end
    hold on
    scatter(length(thresh),thresh(end),'b');
    scatter(length(thresh),av_thresh,'r');
    hold off
end
   
function s = AdaptStairParameters(h,atype)

if isfield(h,'s')
    s = h.s;
else
    s=struct;
    s.out.adaptive = [];
    s.trial = 1;
end

s.p(atype).up=h.Settings.adaptive(atype).updown(1);
s.p(atype).down=h.Settings.adaptive(atype).updown(2);
%s.p(atype).feedback= 1;
s.p(atype).reversals = h.Settings.adaptive(atype).reversals;
s.p(atype).stepsize = h.Settings.adaptive(atype).stepsize;
%s.p(atype).SaveResults =1;
%s.p(atype).tasktype = 1;
%s.p(atype).exppos = 1;
s.p(atype).thresholdtype='Arithmetic';
if isfield(h.Settings.adaptive(atype),'reversalForthresh')
    s.p(atype).reversalForthresh = h.Settings.adaptive(atype).reversalForthresh;
else
    s.p(atype).reversalForthresh = 3;
end
if isfield(h.Settings.adaptive(atype),'steptype')
    s.p(atype).isstep = h.Settings.adaptive(atype).steptype;
else
    s.p(atype).isstep = 0;
end
%s.nafc=3;
%s.NameStepSize='Factor';

% if this is the first run, do some setup
if ~isfield(s,'a') || length(s.a)<atype || isempty(s.a(atype).StimulusLevel)
    setup = 1;
else
    setup=0;
end
if setup
    if length(s.p(atype).reversals) ~= length(s.p(atype).stepsize)
        error('The number of s.reversals and the number of steps must be identical.');
    end
    % here I set the plan of the threshold tracking, i.e., a matrix (two
    % columns) that contains on the left the progressive number of reversals and
    % on the right the corresponding step size
    s.a(atype).expplan = [(1:sum(s.p(atype).reversals))', zeros(sum(s.p(atype).reversals), 1)];
    i=1;
    for j=1:length(s.p(atype).reversals)
        for k=1:s.p(atype).reversals(j)
            s.a(atype).expplan(i, 2)=s.p(atype).stepsize(j);
            i=i+1;
        end
    end
    % here I define the variable 'row of output' that contains all the output values of the
    % function. In the current function the output is updated at the end of the while loop
    % this are the values and the labels
    s.a(atype).rowofoutput = zeros(1, 6);
    s.a(atype).expthresholds = zeros(1, 1);

    %clc
    %input('Press return to begin ', 's');
    %pause(1)
    % indexes for the while loop
    s.a(atype).count_of_n_of_reversals = 0;
    s.a(atype).blockthresholds = zeros(length(s.p(atype).reversalForthresh), 1);
    s.a(atype).n_threshold = 1;
    % variable for the up-s.down
    s.a(atype).n_down = 0;
    % variable for count the positive and negative answers
    s.a(atype).pos = 0;
    s.a(atype).neg = 0;
    s.a(atype).trend = 30;
    s.a(atype).StimulusLevel = h.Settings.adaptive(atype).startinglevel;
    s.a(atype).actualstep = s.a(atype).expplan(1, 2);
end