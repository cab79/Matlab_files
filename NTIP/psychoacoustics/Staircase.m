function [output, expthresholds] = Staircase()

% This is the main function

clear all
clc

tracking = 'Staircase';
s = StaircaseParameters;    %obtaining the staircase procedure parameters

if ~isempty(s)
    if length(s.reversals) ~= length(s.stepsize)
        error('The number of s.reversals and the number of steps must be identical.');
    end
    % here I set the plan of the threshold tracking, i.e., a matrix (two
    % columns) that contains on the left the progressive number of s.reversals and
    % on the right the corresponding step size
    expplan = [(1:sum(s.reversals))', zeros(sum(s.reversals), 1)];
    i=1;
    for j=1:length(s.reversals)
        for k=1:s.reversals(j)
            expplan(i, 2)=s.stepsize(j);
            i=i+1;
        end
    end
    % here I define the variable rowof output that contains all the output values of the
    % function. In the current function the output is updated at the end of the while loop
    % this are the values and the labels
    output = [];
    rowofoutput = zeros(1, 6);
    expthresholds = zeros(s.nblocks, 1);
    % BEGINNING OF THE CICLE WHILE FOR CALCULATING A SINGLE THRESHOLD VALUE
    for block=1:s.nblocks
        clc
        input('Press return to begin the block ', 's');
        clc
        fprintf ('Block number %1.0f\n\n', block);
        pause(1)
        % indexes for the while loop
        count_of_n_of_reversals = 0;
        trial = 1;
        blockthresholds = zeros(length(s.reversalForthresh), 1);
        n_threshold = 1;
        % variable for the up-s.down
        n_down = 0;
        % variable for count the positive and negative answers
        pos = 0;
        neg = 0;
        trend = 30;
        StimulusLevel = s.startinglevel;
        actualstep = expplan(1, 2);
        while count_of_n_of_reversals < sum(s.reversals);
            pause(0.5)
            fprintf('[%1.0f] ', trial);
            % here I get the answer from the simulated listener
            fun = [s.experiment,'(' num2str(s.standard) ', StimulusLevel ,' num2str( s.nafc)   ')'];
            [CorrectAnswer, Question] =  eval(fun);
            s.SubjectAccuracy(trial)= EvaluateAnswer(CorrectAnswer,s.feedback,Question);   %evaluing the subject answer (right or wrong)
            % UPDATE THE ROWOFOUTPUT
            rowofoutput (1, 1) = block;
            rowofoutput (1, 2) = trial;
            rowofoutput (1, 3) = StimulusLevel;
            rowofoutput (1, 4) = s.SubjectAccuracy(trial);
            % here I upddate the count for the up-s.down motion
            if s.SubjectAccuracy(trial) == 1
                n_down = n_down + 1;
                if n_down == s.down
                    n_down = 0;
                    pos = 1;
                    trend = 1;
                    % here I update the count of the number of s.reversals and
                    % corresponding stepsize
                    if pos ==1 && neg == -1
                        count_of_n_of_reversals = count_of_n_of_reversals + 1;
                        % calculate the threshold
                        blockthresholds(n_threshold)=(StimulusLevel + rowofoutput(1, 3))/2;
                        n_threshold = n_threshold + 1;
                        actualstep = expplan(count_of_n_of_reversals, 2);
                        pos = trend;
                        neg = trend;
                    end
                    if s.isstep == 1
                        StimulusLevel = StimulusLevel - actualstep;
                    else
                        StimulusLevel = StimulusLevel / actualstep;
                    end
                end
            else
                neg = -1;
                trend = -1;
                n_down = 0;
                % here I update the count of the number of s.reversals and
                % corresponding stepsize
                if pos ==1 && neg == -1
                    count_of_n_of_reversals = count_of_n_of_reversals + 1;
                    % calculate the threshold
                    blockthresholds(n_threshold)=(StimulusLevel + rowofoutput(1, 3))/2;
                    n_threshold = n_threshold + 1;
                    actualstep = expplan(count_of_n_of_reversals, 2);
                    pos = trend;
                    neg = trend;
                end
                if s.isstep == 1
                    StimulusLevel = StimulusLevel + actualstep;
                else
                    StimulusLevel = StimulusLevel * actualstep;
                end
            end
            % UPDATE THE ROWOFOUTPUT
            rowofoutput (1, 5) = count_of_n_of_reversals;
            rowofoutput (1, 6) = actualstep;
            % update the number of trials
            trial = trial + 1;
            % UPDATE THE GLOBAL OUTPUT VARIABLE
            output = [output; rowofoutput];
        end
        % here I calculate the threshol for the block
        switch s.thresholdtype
            case 'Arithmetic'
                expthresholds(block)=mean(blockthresholds(end-(s.reversalForthresh-1):end));
            case 'Geometric'
                expthresholds(block)=prod(blockthresholds(end-(s.reversalForthresh-1):end))^(1/length(blockthresholds(end-(s.reversalForthresh-1):end)));
            case 'Median'
                expthresholds(block)=median(blockthresholds(end-(s.reversalForthresh-1):end));
            otherwise
                disp('Unknown calculation type.')
        end
        fprintf('Threshold equal to %1.3f\n', expthresholds(block));
        fprintf('Press return to continue\n');
        pause
        fprintf ('\nBLOCK ENDED\n');
        pause(2)
    end
    clc
    fprintf ('\nEXPERIMENT ENDED\n\n');
    s.MATSAVEDATA = output;
    if s.SaveResults
        WriteDataFile(tracking,s.fileout,output,s.nsub,s.name,s.gender,s.age,s.note,expthresholds);
    end
    clear all
end