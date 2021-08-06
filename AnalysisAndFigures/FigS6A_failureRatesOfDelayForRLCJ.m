%%% This script makes the figures showing idiosyncracies for false start
%%% and delay drift failures for the controls
%%%
%%% Adam Smoulder, 4/12/21

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

% Get number of rewards; assume reward labels are same across subjects
rewards = unique([behavior_bySubject{1}.reward]);
nrewards = length(rewards);

%% Get trials by direction, reward, and trial status conditions
getNCondFromBehaviorForMainFigures

%% We'll first make plots for delay drifts (-12) and false starts (-11)
%
% We'll first get the counts and totals for each status for each subject
% separately (different epochs have diff #s of viable reaches)
statusesToCheck = {[1 -22 -23 -34], -11, -12}; % both successes and reach failures are "successful" at getting through the reach epoch
nstatToCheck = length(statusesToCheck);
statusNames = {'Delay Success','False Start','Delay Drift'}; 
nstat_bySubject_byStatus_byReward = nan(nsubjects,nstatToCheck,nrewards); % how often does this status show up?
nfull_bySubject_byStatus_byReward = nan(nsubjects,nstatToCheck,nrewards); % how many trials out of?
for f = 1:nsubjects
    n_byCond = squeeze(sum(n_byCond_bySubject{f}));
    for s = 1:nstatToCheck
        nstat_bySubject_byStatus_byReward(f,s,:) = sum(n_byCond(:,find(ismember(statuses,statusesToCheck{s}))),2); 
        goodStatuses = (statuses <= 10*ceil(statusesToCheck{s}(1)/10)) | (statuses==1); % we use 10s to split up different epochs of failure; anything that made it to this epoch or was success is valid 
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

rewNames = {'Small','Medium','Large','RL/CJ','Jackpot'};
statusColors = [105 190 70; 245 130 48 ; [245 130 48]*0.8 ]/255;
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
    yticks([80 100])
    if f == 1
        ylabel('% of Trials')
    end
    axis([rewards(1)-1 rewards(end)+1 75 100])
    if f ==1
        title('Rare-Large Sessions')
    else
        title('Common-Jackpot Sessions')
    end
end; clear f
set(gcf,'position',[33 670 527 287])
pvals_byRewardXReward_bySubject_byStatus

% Save it!
figname = 'FigS6A_failureRatesOfDelayForRLCJ';
saveas(gcf,['MATLABFigs\' figname])
saveas(gcf,['SVGs\' figname '.svg'])

