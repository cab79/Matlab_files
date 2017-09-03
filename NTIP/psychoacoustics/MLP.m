function MLP()


clear all
clc

tracking = 'MLP';
s = MLParameters;    %obtaining the ML procedure parameters

count=0;
if ~isempty(s)
    for i=1:s.nblocks   %running the blocks
        clc
        input('Press return to begin the block ', 's');
        clc
        fprintf ('Block number %1.0f\n\n', i);
        pause(1)
        for j = 1:s.ntrials %running the trials
            pause(0.5)
            count = count+1;
            if j>1 && ~s.nafc && rand<s.catchtrial
                s.StimulusLevel(j)=s.standard;
            end
            fprintf('[%1.0f] ', j);
            fun = [s.experiment,'(' num2str(s.standard) ', s.StimulusLevel(j) ,' num2str( s.nafc)   ')'];
            [CorrectAnswer, Question] =  eval(fun);
            s.SubjectAccuracy (j)= EvaluateAnswer(CorrectAnswer,s.feedback,Question);   %evaluing the subject answer (right or wrong)
            while j == 1 && s.SubjectAccuracy (j) == 0 && s.repeatft
                fprintf('[%1.0f] ', j);
                fun = [s.experiment,'(' num2str(s.standard) ', s.StimulusLevel(j) ,' num2str( s.nafc)   ')'];
                [CorrectAnswer, Question] =  eval(fun);
                s.SubjectAccuracy (j)= EvaluateAnswer(CorrectAnswer,s.feedback,Question);   %evaluing the subject answer (right or wrong)
            end
            %estimating the stimuli intensity for the successive trial
            [s.TemporaryThreshold(j),s.FA(j)] = FindThreshold (s.p_target, s.StimulusLevel(1:j), ...
                s.SubjectAccuracy (1:j), s.midpoints, s.beta, s.gamma, s.lambda);
            if j < s.ntrials
                s.StimulusLevel(j+1)=s.TemporaryThreshold(j);
            elseif j == s.ntrials
                fprintf('Threshold equal to %1.2f\n', s.TemporaryThreshold(j));
                fprintf('Press return to continue\n');
                pause
            end
            s.MATSAVEDATA(count,:) = [i,j,s.StimulusLevel(j),s.FA(j),s.SubjectAccuracy(j),s.TemporaryThreshold(j)];
        end %end trials
        fprintf ('\nBLOCK ENDED\n');
        pause(2)
    end  %end blocks
    clc
    fprintf ('\nEXPERIMENT ENDED\n\n');
    if s.SaveResults
        WriteDataFile(tracking,s.fileout,s.MATSAVEDATA,s.nsub,s.name,s.gender,s.age,s.note);
    end
    clear all
else
    clear all
end