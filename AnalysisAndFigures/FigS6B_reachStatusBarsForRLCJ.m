%%% This script makes the figures showing overshoot/undershoot rates for
%%% the control experiments.
%%%
%%% Adam Smoulder, 8/21/20 (edit 5/4/21)

% Load subject data behavior files
subjectNames = {'EarlRL','EarlCJ'};
behaviorFnames = {...
    'EarlRL_behaviorFeatures_20210107',...
    'EarlCJ_behaviorFeatures_20210107',...
    };
nsubjects = length(subjectNames);
behavior_bySubject = cell(nsubjects,1);
for f = 1:nsubjects
    load(behaviorFnames{f},'behavior');
    behavior_bySubject{f} = behavior;
    disp(['Loaded subject ' subjectNames{f}])
end; clear f
clear behavior

%% Get trials by direction, reward, and trial status conditions
getNCondFromBehaviorForMainFigures

%% 
% We'll first get the counts and totals for each status for each subject
% separately (different epochs have diff #s of viable reaches)
statusesToCheck = {[1 -34], -22, -23}; % both successes and reach failures are "successful" at getting through the reach epoch
nstatToCheck = length(statusesToCheck);
statusNames = {'Reach Success','Overshoot','Undershoot'}; 
nstat_bySubject_byStatus_byReward = nan(nsubjects,nstatToCheck,nrewards); % how often does this status show up?
nfull_bySubject_byStatus_byReward = nan(nsubjects,nstatToCheck,nrewards); % how many trials out of?
for f = 1:nsubjects
    n_byCond = squeeze(sum(n_byCond_bySubject{f}));
    for s = 1:nstatToCheck
        nstat_bySubject_byStatus_byReward(f,s,:) = sum(n_byCond(:,find(ismember(statuses,statusesToCheck{s}))),2); 
        goodStatuses = ismember(statuses,[statusesToCheck{:}]); % we use 10s to split up different epochs of failure; anything that made it to this epoch or was success is valid 
        nfull_bySubject_byStatus_byReward(f,s,:) = sum(n_byCond(:,goodStatuses),2); % all trials that made it to the relevant epoch are used here
    end; clear s
end; clear f

statRates_bySubject_byStatus_byReward = ...
    nstat_bySubject_byStatus_byReward./nfull_bySubject_byStatus_byReward*100;


%% We can do binomial proportion test and find SEMs for these rates
pvals_byRewardXReward_bySubject_byStatus = nan([nrewards,nrewards,nsubjects,nstatToCheck]);
statRates_mean_bySubject_byStatus_byReward = nan(size(statRates_bySubject_byStatus_byReward));
statRates_sem_bySubject_byStatus_byReward = nan(size(statRates_bySubject_byStatus_byReward));
nboots = 10000;
for f = 1:nsubjects
    for s = 1:nstatToCheck
        eventCount = squeeze(nstat_bySubject_byStatus_byReward(f,s,:));
        totalCount = squeeze(nfull_bySubject_byStatus_byReward(f,s,:));
        pvals_byRewardXReward_bySubject_byStatus(:,:,f,s) = ...
            binomialProportionTest(eventCount,totalCount);
        [~,bootMean,~,~,~,~,bootSem] = ...
            bootstrapBinaryEvent(eventCount,totalCount,nboots);
        statRates_mean_bySubject_byStatus_byReward(f,s,:) = bootMean./totalCount'*100;
        statRates_sem_bySubject_byStatus_byReward(f,s,:) = bootSem./totalCount'*100;
    end; clear s
    disp(['Completed bootstrapping for ' subjectNames{f}])
end; clear f



%% Show plots for each subject

rewNames = {'Small','Medium','Large','CJ/RL','Jackpot'};
statusColors = [105 190 70; 200 200 200 ; 50 50 50 ]/255;
lw = 0.5;
ms = 10;
fs = 10;

figure
for f = 1:nsubjects
    subplot(1,nsubjects,f); hold on
    curSuccessFOR = squeeze(statRates_bySubject_byStatus_byReward(f,1,:)); % FOR = fraction of reaches
    curFSFOR = squeeze(statRates_bySubject_byStatus_byReward(f,2,:));
    curDDFOR = squeeze(statRates_bySubject_byStatus_byReward(f,3,:));
    
    p = bar([curSuccessFOR curDDFOR curFSFOR ],'stacked');
    p(1).FaceColor = statusColors(1,:);
    p(2).FaceColor = statusColors(2,:);
    p(3).FaceColor = statusColors(3,:);
    
    curSuccessSEM = squeeze(statRates_sem_bySubject_byStatus_byReward(f,1,:));
%     errorbar(rewards,curSuccessFOR,curSuccessSEM,'k.-','linewidth',lw,'markersize',ms)
    plot(rewards,curSuccessFOR,'k.-','linewidth',lw,'markersize',ms)
    
    set(gca,'fontname','arial')
    set(gca,'fontsize',fs)
    xticks(1:nrewards)
    xticklabels(rewNames)
    xtickangle(45)
    yticks([50 75 100])
    if f == 1
        ylabel('% of Trials')
    end
    axis([rewards(1)-1 rewards(end)+1 50 100])
    if f ==1
        title('Rare-Large Sessions')
    else
        title('Common-Jackpot Sessions')
    end
end; clear f
set(gcf,'position',[33 670 527 287])

pvals_byRewardXReward_bySubject_byStatus

% Save it!
figname = 'FigS6B_reachStatusBarsForRLCJ';
saveas(gcf,['MATLABFigs\' figname])
saveas(gcf,['SVGs\' figname '.svg'])

% Fascinatingly, it appears Earl's reach choke modality is now
% overshoots??? Man, good thing we saved this exp for last, CJ really did
% mess with him...
