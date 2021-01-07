%% This script wraps the count-by-cond-from-behavior function. We use 
%  this to consolidate some of the failure methods based on the analyses
%  that we're doing. This script specifically is for the figures of the
%  main paper, where:
%  - quitouts and misstarts (no real trial) are excluded ([-20 -13 -14 0] = 0)
%  - No attempts, wild, and early returns are also removed due to infrequency ([-21 -26 -33] = 0)
%  - the few -111s we have are misstarts (-111 = 0)
%  - we'll remove any excluded trials (0s removed)
%  - overshoots = targ hold overshoots = scuff ([-22 -32 -31] = -22)
%  - undershoots = slow ([-23 -24 -25] = -23)
%  - target hold failures = jitter also ([-34 -35] = -34)
%
% Assumes that behavior_bySubject is in the workspace, along with nsubjects
% and subjectNames
%
% Adam Smoulder, 7/24/20 (edit 10/13/20)

% If you want to include the excluded statuses in the main analyses (e.g.
% like Fig 1B), just put their #s into the "statuses" variable below

n_byCond_bySubject = cell(nsubjects,1);
dirLabels_bySubject = cell(nsubjects,1);
rewardLabels_bySubject = cell(nsubjects,1);
trialStatusLabels_bySubject = cell(nsubjects,1);

statuses = [1 -11 -12 -22 -23 -34];
nstatuses = length(statuses);
for f = 1:nsubjects % for each subject
    % First, run the main count on the behavior
    [n_byCondOrig, dirLabels, rewardLabels, trialStatusLabels, dirs,...
    rewards, curStatuses] = getNByCondFromBehavior(behavior_bySubject{f});
    ndirs = length(dirs);
    nrewards = length(rewards);

    % Then, we do all the status combining and such that's desired (see
    % descriptions above). We also skip all the 0 (or other skipped things)
    n_byCond = zeros(ndirs, nrewards, nstatuses);
    for s = 1:nstatuses
        if statuses(s)==-22 % overshoots
            n_byCond(:,:,s) = sum(n_byCondOrig(:,:,ismember(curStatuses,[-22 -31 -32])),3);
        elseif statuses(s)==-23 % undershoots
            n_byCond(:,:,s) = sum(n_byCondOrig(:,:,ismember(curStatuses,[-23 -24 -25])),3);
        elseif statuses(s)==-34 % target hold drifts
            n_byCond(:,:,s) = sum(n_byCondOrig(:,:,ismember(curStatuses,[-34 -35])),3);
        elseif ~ismember(statuses(s),curStatuses) % the status we're looking for isn't present in this data; skip
            continue
        else % anything else
            n_byCond(:,:,s) = n_byCondOrig(:,:,find(curStatuses==statuses(s),1));
        end
    end; clear s
    
    % Save this to the overall cell array, along with the labels
    n_byCond_bySubject{f} = n_byCond;
    dirLabels_bySubject{f} = dirLabels;
    rewardLabels_bySubject{f} = rewardLabels;
    trialStatusLabels_bySubject{f} = trialStatusLabels;
    
    disp(['Completed getting condition counts for ' subjectNames{f}])
end; clear f n_byCond dirLabels rewardLabels trialStatusLabels
