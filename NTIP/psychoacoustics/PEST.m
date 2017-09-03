function [threshold, output] = PEST()

% This is the main function

clear all
% clc

tracking = 'PEST';
s = PestParameters;    %obtaining the ML procedure parameters

if ~isempty(s)
    % the variables to store resultsss
    rowofoutput = zeros(1, 5);
    output = [];
    threshold = zeros(s.nblocks, 1);
    for block=1:s.nblocks
        clc
        input('Press return to begin the block ', 's');
        clc
        fprintf ('Block number %1.0f\n\n', block);
        pause(1)
        StimulusLevel = s.startinglevel;
        step=s.startingstepsize;
        trial=0;
        historycorrect = '';
        doubled = 0;
        while step > s.finalstepsize
            n_correct = 0;
            runtrial=0;
            % LIKELIHOOD TEST
            while (s.p_target * runtrial) - s.waldfactor < n_correct && n_correct < (s.p_target * runtrial) + s.waldfactor
                pause(0.5)
                trial = trial + 1;
                runtrial = runtrial + 1;
                fprintf('[%1.0f] ', trial);
                fun = [s.experiment,'(' num2str(s.standard) ', StimulusLevel ,' num2str( s.nafc)   ')'];
                [CorrectAnswer, Question] =  eval(fun);
                s.SubjectAccuracy(trial)= EvaluateAnswer(CorrectAnswer,s.feedback,Question);   %evaluing the subject answer (right or wrong)
                n_correct = n_correct + s.SubjectAccuracy(trial);
                rowofoutput(1)=block;
                rowofoutput(2)=trial;
                rowofoutput(3)=StimulusLevel;
                rowofoutput(4)=s.SubjectAccuracy(trial);
                rowofoutput(5)=step;
                output = [output; rowofoutput];
            end
            
            if n_correct >= ((runtrial*s.p_target)+s.waldfactor)
                historycorrect = [historycorrect, 'd'];
            elseif n_correct <= ((runtrial*s.p_target)-s.waldfactor)
                historycorrect = [historycorrect, 'u'];
            end
            
            % PEST: RULE N. 1
            if length(historycorrect) > 1
                if historycorrect(end-1)~= historycorrect(end)
                    step=step/2;
                end
            end
            
            % RULE N. 3 AND 4
            if length(historycorrect) > 2
                if length(historycorrect) == 3 %3 consecutive steps in the same direction, if they occur on the first steps, will not be
                    %preceeded by a reversal, nonetheless, I think the step size should be doubled in this case
                    if strcmp(historycorrect(end-2:end),'ddd')|| strcmp(historycorrect(end-2:end),'uuu')
                        step = step * 2;
                        doubled=1;
                    end
                end
                
                if length(historycorrect) > 3
                    if strcmp(historycorrect(end-3:end),'uddd')|| strcmp(historycorrect(end-3:end),'duuu')
                        if doubled
                            doubled=0;
                        else
                            step = step * 2;
                            doubled=1;
                        end
                    end
                    if strcmp(historycorrect(end-3:end), 'uuuu') || strcmp(historycorrect(end-3:end), 'dddd')
                        step = step * 2;
                        doubled=1;
                    end
                end
            end
            
            % EXCEPTION TO RULE N. 3 (TAYLOR & CREELMAN, 1967)
            if step > s.maxstepsize
                step = s.maxstepsize;
            end
            
            if strcmp(historycorrect(end),'u')
                StimulusLevel = StimulusLevel + step;
            else
                StimulusLevel = StimulusLevel - step;
            end
            
            if StimulusLevel <= s.minlevel
                StimulusLevel = s.minlevel + s.finalstepsize;
            end
        end
        threshold(block) = StimulusLevel;
        fprintf('Threshold equal to %4.3f\n', threshold(block));
        fprintf('Press return to continue\n');
        pause
        fprintf ('\nBLOCK ENDED\n');
        pause(2)
    end
    clc
    fprintf ('\nEXPERIMENT ENDED\n\n');
    s.MATSAVEDATA = output;
    if s.SaveResults
        WriteDataFile(tracking,s.fileout,output,s.nsub,s.name,s.gender,s.age,s.note,threshold);
    end
    clear all
end