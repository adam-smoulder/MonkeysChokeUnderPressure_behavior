%%% This script makes the figures showing idiosyncracies for target hold
%%% failures
%%%
%%% Adam Smoulder, 8/21/20 (edit 10/13/20)

% Load subject data behavior files
loadBehaviorDataForMainFigures

%% Get trials by direction, reward, and trial status conditions
getNCondFromBehaviorForMainFigures

%% Make a plot for the successes and target hold failures
%
% We'll first get the counts and totals for each status for each subject
% separately (different epochs have diff #s of viable reaches)
statusesToCheck = {1,-34}; % both successes and reach failures are "successful" at getting through the reach epoch
nstatToCheck = length(statusesToCheck);
statusNames = {'Success','Drift'};
nstat_bySubject_byStatus_byReward = nan(nsubjects,nstatToCheck,nrewards); % how often does this status show up?
nfull_bySubject_byStatus_byReward = nan(nsubjects,nstatToCheck,nrewards); % how many trials out of?
for f = 1:nsubjects
    n_byCond = squeeze(sum(n_byCond_bySubject{f}));
    for s = 1:nstatToCheck
        nstat_bySubject_byStatus_byReward(f,s,:) = sum(n_byCond(:,find(ismember(statuses,statusesToCheck{s}))),2); 
        goodStatuses = ismember(statuses,[1 -34]); % we use 10s to split up different epochs of failure; anything that made it to this epoch or was success is valid 
        nfull_bySubject_byStatus_byReward(f,s,:) = sum(n_byCond(:,goodStatuses),2); % all trials that made it to the relevant epoch are used here
    end; clear s
end; clear f

statRates_bySubject_byStatus_byReward = ...
    nstat_bySubject_byStatus_byReward./nfull_bySubject_byStatus_byReward*100;

% Quick aside - we need to know the overall % of hold failures that occurs
% out of successes for the caption of the endpoints figure
chanceOfFailureDuringTargetHoldForN = ...
    sum(nstat_bySubject_byStatus_byReward(1,2,:))/sum(nfull_bySubject_byStatus_byReward(1,2,:))*100


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

rewNames = {'Small','Medium','Large','Jackpot'};
statusColors = [105 190 70; [255 100 255] ]/255;
lw = 0.5;
ms = 10;
fs = 10;

figure
for f = 1:nsubjects
    subplot(1,nsubjects,f); hold on
    curSuccessFOR = squeeze(statRates_bySubject_byStatus_byReward(f,1,:)); % FOR = fraction of reaches
    curDriftFOR = squeeze(statRates_bySubject_byStatus_byReward(f,2,:));
%     curJitterFOR = squeeze(statRates_bySubject_byStatus_byReward(f,3,:));
    
    p = bar([curSuccessFOR curDriftFOR ],'stacked');
    p(1).FaceColor = statusColors(1,:);
    p(2).FaceColor = statusColors(2,:);
%     p(3).FaceColor = statusColors(3,:);
    
    curSuccessSEM = squeeze(statRates_sem_bySubject_byStatus_byReward(f,1,:));
%     errorbar(rewards,curSuccessFOR,curSuccessSEM,'k.-','linewidth',lw,'markersize',ms)
    plot(rewards,curSuccessFOR,'k.-','linewidth',lw,'markersize',ms)
    
    set(gca,'fontname','arial')
    set(gca,'fontsize',fs)
    xticks(1:nrewards)
    xticklabels(rewNames)
    xtickangle(45)
    yticks([85 100])
    if f == 1
        ylabel('% of Trials')
    end
    axis([rewards(1)-1 rewards(end)+1 81 100])
    title(['Monkey ' subjectNames{f}(1)])
end; clear f
set(gcf,'position',[33 670 790 287])

% Save it!
figname = 'FigS5C_failureRatesOfTargetHold';
saveas(gcf,['MATLABFigs\' figname])
saveas(gcf,['SVGs\' figname '.svg'])

