%%% This script gets all of the values for epoch success rates and specific
%%% failure modes.
%%%
%%% Adam Smoulder, 10/13/20

% Load subject data behavior files
loadBehaviorDataForMainFigures

%% Get trials by direction, reward, and trial status conditions
getNCondFromBehaviorForMainFigures

%% Epoch success rates first
n_bySubject_byReward_byStatus = nan(nsubjects,nrewards,nstatuses);
for f = 1:nsubjects
    n_bySubject_byReward_byStatus(f,:,:) = sum(n_byCond_bySubject{f}); % marginalizing over reach directions
end

epochSuccStatuses = {[1 -22 -23 -33],[1 -34],[1]};
epochFailStatuses = {[-11 -12],[-22 -23],[-34]};
nepochs = length(epochSuccStatuses);
n_bySubject_byReward_byEpoch = nan(nsubjects,nrewards,nepochs);
nsucc_bySubject_byReward_byEpoch = nan(nsubjects,nrewards,nepochs);
for e = 1:nepochs
    nsucc_bySubject_byReward_byEpoch(:,:,e) = sum(n_bySubject_byReward_byStatus(:,:,ismember(statuses,epochSuccStatuses{e})),3);
    n_bySubject_byReward_byEpoch(:,:,e) = sum(n_bySubject_byReward_byStatus(:,:,ismember(statuses,epochFailStatuses{e})),3)...
                                          + nsucc_bySubject_byReward_byEpoch(:,:,e);
end; clear e

succRate_bySubject_byReward_byEpoch = ...
    100*nsucc_bySubject_byReward_byEpoch./n_bySubject_byReward_byEpoch;

%% Get p-values
pvals_byRewardXReward_bySubject_byEpoch = nan(nrewards,nrewards,nsubjects,nepochs);
for f = 1:nsubjects
    for e = 1:nepochs
        eventCounts = nsucc_bySubject_byReward_byEpoch(f,:,e);
        totalCounts = n_bySubject_byReward_byEpoch(f,:,e);
        pvals_byRewardXReward_bySubject_byEpoch(:,:,f,e) = binomialProportionTest(eventCounts,totalCounts);
    end; clear e
end; clear f

%% For display, this shows each subject's S->L and L->J change for each epoch
diff(permute(succRate_bySubject_byReward_byEpoch(:,[1 3 4],:),[3 2 1]),[],2)

%% For display, this shows each p-value individually
for f = 1:nsubjects
    disp(['Subject ' subjectNames{f}(1)])
    for e = 1:nepochs
        disp(['Epoch ' num2str(e) ', S->L p = ' num2str(pvals_byRewardXReward_bySubject_byEpoch(1,3,f,e))])
        disp(['Epoch ' num2str(e) ', L->J p = ' num2str(pvals_byRewardXReward_bySubject_byEpoch(3,4,f,e))])
    end; clear e
end; clear f







%% Now individual failure rates
failEpochs = [1 1 2 2 3]; % epochs for different failures
failRate_bySubject_byReward_byStatus = 100*...
    n_bySubject_byReward_byStatus(:,:,2:end)./n_bySubject_byReward_byEpoch(:,:,failEpochs);

%% Get p-values
pvals_byRewardXReward_bySubject_byStatus = nan(nrewards,nrewards,nsubjects,nepochs);
for f = 1:nsubjects
    for a = 1:nstatuses-1
        eventCounts = n_bySubject_byReward_byStatus(f,:,a+1);
        totalCounts = n_bySubject_byReward_byEpoch(f,:,failEpochs(a));
        pvals_byRewardXReward_bySubject_byStatus(:,:,f,a) = binomialProportionTest(eventCounts,totalCounts);
    end; clear e
end; clear f


%% For display, this shows each subject's S->L and L->J change for each epoch
diff(permute(failRate_bySubject_byReward_byStatus(:,[1 3 4],:),[3 2 1]),[],2)

%% For display, this shows each p-value individually
for f = 1:nsubjects
    disp(['Subject ' subjectNames{f}(1)])
    for a = 1:nstatuses-1
        disp(['Status ' num2str(statuses(a+1)) ', S->L p = ' num2str(pvals_byRewardXReward_bySubject_byStatus(1,3,f,a))])
        disp(['Status ' num2str(statuses(a+1)) ', L->J p = ' num2str(pvals_byRewardXReward_bySubject_byStatus(3,4,f,a))])
    end; clear e
end; clear f
