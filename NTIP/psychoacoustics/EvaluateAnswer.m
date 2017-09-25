function subject_answer = EvaluateAnswer(ra,feedback,q)
%
% subject_answer = EvaluateAnswer(ra, q)
% 
% This function echoes on screen the question the subject has to answer.
% In yes/no tasks it codes the subject's answer into 0 (i.e., no) or 1
% (i.e., yes). In nafc tasks it also matches the subject's response with
% the actual correct response.
% - RA: the answer given by the subject
% - Q: a text with the question the subjects has to answer

if nargin<3
    q='';
end;


disp(q)
user_input= str2double(getkey (1,'non-ascii'));


if ra == user_input
    subject_answer = 1;
    if feedback
        disp('CORRECT');
    end;
else
    subject_answer = 0;
    if feedback
        disp('WRONG');
    end;
end;