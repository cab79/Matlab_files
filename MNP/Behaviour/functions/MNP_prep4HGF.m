function [D,S] = MNP_prep4HGF(D,S,fitsim)

for d = 1:length(D)
    % associative learning: D.HGF.u indicates outcome and cues
    sig=D(d).Sequence.signal(1:2,:); 
    D(d).HGF.u(1,sig(2,:)==1) = 0; % outcome
    D(d).HGF.u(1,sig(2,:)==2) = 1; % outcome
    D(d).HGF.u(2,:) = 1;%sig(2,:); % outcome types
    D(d).HGF.u(3,:) = sig(1,:); % cues
    D(d).HGF.u=D(d).HGF.u';

    % if not simulating, gather response data
    if S.fitsim==1

        % create response matrix
        D(d).HGF.y=nan(size(D(d).Sequence.signal,2),2);

        %choices
        if any(strcmp(S.resp_modelspec.responses,'Ch'))
            D(d).HGF.y(:,1) = D(d).Processed.presssignal; % BINARY response
            D(d).HGF.y(D(d).HGF.y(:,1)==1)=0;
            D(d).HGF.y(D(d).HGF.y(:,1)==2)=1;
        end

        % response times
        if any(strcmp(S.resp_modelspec.responses,'RT'))
            RT=D(d).Output.RT;
            RT(RT<0.5)=nan; % don't consider RTs less than 500ms
            D(d).HGF.y(D(d).Output.presstrial,2) = log(RT);
        end

        %S.use_y_col = find([any(strcmp(S.resp_modelspec.responses,'Ch')), any(strcmp(S.resp_modelspec.responses,'RT')),any(strcmp(S.resp_modelspec.responses,'EEG'))]);
    end
end

%% OLD

    % stimulus - binary (0 and 1)
%     D.HGF.u=D.Sequence.signal(2,:); % second stimulus
%     D.HGF.u(D.HGF.u==1)=0;
%     D.HGF.u(D.HGF.u==2)=1;
%     D.HGF.u(2,:)=1;
%     D.HGF.u=D.HGF.u';
      %if isempty(sim_param)
    %     D.HGF.y = D.Processed.presssignal; % BINARY response
    %     D.HGF.y(D.HGF.y==1)=1;
    %     D.HGF.y(D.HGF.y==2)=0;
    %     D.HGF.y=D.HGF.y';
      %end

%     % associative learning: D.HGF.u indicates pairings
%     sig=D.Sequence.signal(1:2,:); 
%     D.HGF.u(sig(1,:)==sig(2,:)) = 0;
%     D.HGF.u(sig(1,:)~=sig(2,:)) = 1;
%     D.HGF.u(2,:)=1;
%     D.HGF.u=D.HGF.u';
    %if isempty(sim_param)
    %     % associative learning - responses indicate pairings
    %     D.HGF.y=[]
    %     ysig=D.Processed.presssignal;
    %     D.HGF.y(ysig(1,:)==sig(1,:)) = 1;
    %     D.HGF.y(ysig(1,:)~=sig(1,:)) = 0;
    %     D.HGF.y(isnan(ysig))=nan;
    %     D.HGF.y=D.HGF.y';
    %end
    
    