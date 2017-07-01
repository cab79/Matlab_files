function condevents = EEG2condevents(EEG, stimmarker, eventtypes,use_etype)

% from an EEGLAB structure, returns a index of conditions of interest
% (length of trials) according to specified stimulus markers

% currently tailored to CRPS digit perception study - requires generifying

%% initiate
if isempty(stimmarker); stimmarker='STIM';end;

%% determine number of conditions
if length(use_etype)==2
    condtypes = combvec(eventtypes{2,use_etype(2)},eventtypes{2,use_etype(1)});
    condtypes(1:2,:) = condtypes([2 1],:);
else
    condtypes = eventtypes{2,use_etype};
end
no_cond = size(condtypes,2);

events = zeros(1,size(EEG.data,3));
%j=0;
for i = 1:size(EEG.data,3)
    if sum(strcmp(EEG.epoch(i).eventtype(1,:),stimmarker))>0
        ind = find(strcmp(EEG.epoch(i).eventtype(1,:),stimmarker));
        for et = 1:length(use_etype)
            findnum=[];
            if iscell(EEG.epoch(i).eventinit_index)
                findnum = find(strcmp(eventtypes{1, use_etype(et)},EEG.epoch(i).eventcodes{1,ind(1)}(:,1)));
                try 
                    events(et,i) = EEG.epoch(i).eventcodes{1,ind(1)}{findnum(end),2};
                catch
                    i
                    et
                    pause
                end
            else 
                findnum = find(strcmp(eventtypes{1, use_etype(et)},EEG.epoch(i).eventcodes(:,ind(1))));
                events(et,i) = EEG.epoch(i).eventcodes{findnum(end),2};
            end
        end
    end
end

if any(strcmp(eventtypes(1,use_etype),'CNUM'))
    events(strcmp(eventtypes(1,:),'CNUM'),:) = events(strcmp(eventtypes(1,:),'CNUM'),:)>0; % change CNUM events to just differentiate change or no change, rather than graduated change
end

eventuniq = unique(events','rows');

if any(strcmp(eventtypes(1,use_etype),'FNUM'))
    if all(ismember(eventuniq(:,1),condtypes(1,:)+5)) % if Right hand
        condtypes(1,:) = condtypes(1,:)+5;
        condtypes(1,:)=sort(condtypes(1,:),'descend');
    end
end

condevents=[];
for c=1:no_cond
    condeventindx = find(ismember(events',condtypes(:,c)','rows'));
    condevents(1,condeventindx) = c*ones(1,length(condeventindx)); 
end
