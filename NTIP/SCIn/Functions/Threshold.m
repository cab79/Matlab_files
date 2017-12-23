function h = Threshold(h)

if ~isfield(h.Settings,'threshold')
    return
end

if isfield(h,'s')
    s=h.s;
else
    s = ThresholdParameters(h); % see function at the end
end

if ~isfield(s,'lasttrialrun')
    s.lasttrialrun = 0;
end
if ~isfield(s,'lastpresstrial')
    s.lastpresstrial = 0;
end


%% RUN


%if strcmp(opt,'responded')
    
%else
%    s.lasttrialrun = h.i;
%end
    
% evaluate the subject's response
resi=0;
if ~isempty(h.out.presstrial)
    % only run once per trial if pressed
    if ~(s.lastpresstrial == h.out.presstrial(end))
        s.lastpresstrial = h.out.presstrial(end);
    else
        % only run once per trial if NOT pressed
        if ~(s.lasttrialrun == h.i)
            s.lasttrialrun = h.i;
        else
            return
        end
    end
    if h.out.presstrial(end) == h.i
        %presstrial=trls_pressed(end);
        pressbutton = h.out.pressbutton(end);
        resi = find(strcmp(pressbutton(end),h.Settings.buttonopt)); % which button was pressed?
    end
else
    % only run once per trial if NOT pressed
    if ~(s.lasttrialrun == h.i)
        s.lasttrialrun = h.i;
    else
        return
    end
end
disp('Running Threshold')
if resi
%s.SubjectAccuracy(s.adaptive.trial)= EvaluateAnswer(CorrectAnswer,s.feedback,Question);   %evaluing the subject answer (right or wrong)
    s.change= -1;
else
    s.change= 1;
end

% UPDATE THE ROWOFOUTPUT
s.rowofoutput (1, 1) = s.trial;
s.rowofoutput (1, 2) = s.StimulusLevel;
s.rowofoutput (1, 3) = s.actStimulusLevel;
s.rowofoutput (1, 4) = s.change;
% here I upddate the count for the up-s.down motion
if s.change == 1
    if s.isstep == 1
        s.StimulusLevel = s.StimulusLevel + s.actualstep;
    else
        s.StimulusLevel = s.StimulusLevel * s.actualstep;
    end
else
    s.StimulusLevel = h.Settings.threshold.startinglevel;
    %if s.isstep == 1
    %    s.StimulusLevel = s.StimulusLevel + s.actualstep;
    %else
    %    s.StimulusLevel = s.StimulusLevel * s.actualstep;
    %end
end

% additional attenutation from GUI or settings
if isfield(h,'vol_atten')
    try
        inten_atten = str2double(get(h.vol_atten,'string'));
    catch
        inten_atten = str2double(h.vol_atten);
    end
else
    inten_atten = h.Settings.atten; 
end

s.actStimulusLevel = inten_atten+s.StimulusLevel;

% UPDATE THE ROWOFOUTPUT
s.trial = s.trial + 1;

%disp(['length_blockthresh = ' num2str(length(s.blockthresholds))]);
disp(['level = ' num2str(s.actStimulusLevel)]);

% threshold calculation
s.expthresholds=0;
if ~isempty(s.out.threshold)
    presstrialsall = find(s.out.threshold(:,4)==-1);
    if any(presstrialsall)
        s.expthresholds = mean(s.out.threshold(presstrialsall-1,3));
        fprintf('Threshold equal to %1.3f\n', s.expthresholds);
    end
end
s.rowofoutput (1, 5) = s.expthresholds;

% UPDATE THE GLOBAL OUTPUT VARIABLE
s.out.threshold = [s.out.threshold; s.rowofoutput];
h.s =s;
h.out.threshold = s.out.threshold;
%fprintf('Press return to continue\n');
%pause
%fprintf ('\nBLOCK ENDED\n');
%pause(2)
end

function s = ThresholdParameters(h)

s.stepsize = 2;
s.SaveResults =1;
s.numForthresh = 3;
s.isstep = 1;


% if this is the first run, do some setup
if ~isfield(s,'trial')
    s.out.threshold = [];
    s.rowofoutput = zeros(1, 6);
    s.expthresholds = zeros(1, 1);
    s.trial = 1;
    s.blockthresholds = zeros(length(s.numForthresh), 1);
    s.StimulusLevel = h.Settings.threshold.startinglevel;
    s.actStimulusLevel = NaN;
    s.actualstep = h.Settings.threshold.step;
end

end