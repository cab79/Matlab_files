function F = fAnova( Dataset )
%% INPUT: First 3 columns - coordinates, 4th column - group index   
%% OUTPUT F-statistic
[GroupLevels,~,GroupIDXs]=unique(Dataset(:,4)); 
% Get number of different Groups
NumberOfGroups=numel(GroupLevels); 
% Get number of observations, equal to number of rows
NumberOfObservations=size(Dataset,1); 
% Count number of observations in each group. It is the same as nItems, but
% lets imagine we work not with simulated data
NumberOfObservationsInEachGroup=accumarray(GroupIDXs,Dataset(:,4),[],@(x) numel(x)); 
% Calculate pair-wise distances
PairWiseDistances=pdist(Dataset(:,1:3));
% make distance matrix
DistMatrix=squareform(pdist(Dataset(:,1:3)));
% leave only low triangle
DistMatrix=tril(DistMatrix);
% Get grouping matrix. 0 when observations from different groups, ones when
% from same group.
GroupMatrix=repmat(Dataset(:,4)',NumberOfObservations,1)==repmat(Dataset(:,4),1,NumberOfObservations);
% Calculate matrix with number of observations per group. Equal to
% GroupMatrix, but instead of 1 for each group we have number of
% observations.
NumberOfObservationsMatrixRedundant=[];
for iGroup=1:NumberOfGroups
   NumberOfObservationsMatrixRedundant=[NumberOfObservationsMatrixRedundant,repmat(NumberOfObservationsInEachGroup(iGroup),NumberOfObservations,NumberOfObservationsInEachGroup(iGroup))];
end %iGrooup
% Make zeros for cells from different groups
NumberOfObservationsMatrix=NumberOfObservationsMatrixRedundant.*GroupMatrix;
%Redundant
% Get SS_T
SSt=sum(sum(DistMatrix.^2))/NumberOfObservations;
% Get SS_W

SSw=nansum(nansum(((DistMatrix.*GroupMatrix).^2)./NumberOfObservationsMatrix));

% Get SS_A
SSa=SSt-SSw;

% Get F
F=abs((SSa/(NumberOfGroups-1))/(SSw/(NumberOfObservations-NumberOfGroups)));

end % Function