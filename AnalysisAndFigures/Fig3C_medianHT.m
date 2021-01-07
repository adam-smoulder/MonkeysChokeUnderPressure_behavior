%%% This script makes the data figures showing the different timings -
%%% reaction time, halfway time, and homing time. We only show homing time
%%% bc that's what's used in the paper
%%%
%%% Adam Smoulder, 7/28/20 (edit 8/7/20)

% Load subject data behavior files
loadBehaviorDataForMainFigures

% Get number of rewards; assume reward labels are same across subjects
rewards = unique([behavior_bySubject{1}.reward]);
nrewards = length(rewards);

%% Find the timings - set parameters for homing time endpoints below
% It seems like the best parameters to use to maximize included trials
% while staying kosher are 1mm away from the target with up to 150ms after
% fail time. The effects we're looking at honestly only get stronger if
% you're more exclusionary, so we're showing conservative results if
% anything.
homingStopDistance = 1; % mm - how far from end target to stop homing time calculation
homingStopAddedTime = 150; % ms - how many ms beyond fail time to allow for homing time calc

reactionTime_bySubject_byReward = cell(nsubjects,nrewards);
midreachTime_bySubject_byReward = cell(nsubjects,nrewards);
homingTime_bySubject_byReward = cell(nsubjects,nrewards);
reachTime_bySubject_byReward = cell(nsubjects,nrewards);

for f = 1:nsubjects
    behavior = behavior_bySubject{f};
    startTargRad = behavior(1).task_parameters.startTargRadius;
    reachTargRad = behavior(1).task_parameters.endTargRadius;
    failTime = behavior(1).task_parameters.failTime;
    targetDistance = behavior(1).task_parameters.targetDistance;
    homingStartDistance = (targetDistance-startTargRad-reachTargRad)*1/3+reachTargRad % use 1/3 left to target
    homingEndDistance = homingStopDistance+reachTargRad;
    maxHomingEndTime = failTime+homingStopAddedTime;
    
    % Remove delay failures, quit outs, no-attempts
    goodBeh = behavior(~ismember([behavior.trialStatusLabels],[0 -111 -11 -12 -13 -14 -20 -21 -26]));
    ngood = length(goodBeh);
    
    % Get all the timings
    reactionTime = nan(ngood,1);
    midreachTime = nan(ngood,1);
    homingTime = nan(ngood,1);
    for i = 1:ngood
        % Get reaction time
        goCueTime = goodBeh(i).t3_goCueTimes;
        reactionTime(i) = goodBeh(i).t4c_reactionTimes_exit-goCueTime;
        
        % Get halfway time
        kinTime = goodBeh(i).kinematics_updated.time;
        distsFromEnd = goodBeh(i).kinematics_updated.distanceFromEndTarget;
        upperMidInd = find(distsFromEnd <= homingStartDistance,1);
        lowerMidInd = upperMidInd-1;
        if isempty(upperMidInd)
            disp(['No mid-reach completed, subject ' subjectNames{f}(1) ' goodBeh ind ' num2str(i)])
            continue
        else
            midreachTime(i) = interp1(distsFromEnd(lowerMidInd:upperMidInd),...
                kinTime(lowerMidInd:upperMidInd),homingStartDistance)...
                -reactionTime(i)-goCueTime;
        end
        
        % If halfway time is beyond trial limits, exclude
        if midreachTime(i)+reactionTime(i) > failTime
            disp(['Super late mid-reach completion, subject ' subjectNames{f}(1) ' goodBeh ind ' num2str(i), ', HT+RT =  ' num2str(midreachTime(i)+reactionTime(i))])
            midreachTime(i) = nan;
            continue
        end
        
        % Get homing time
        upperHomingInd = find(distsFromEnd <= homingEndDistance,1);
        lowerHomingInd = upperHomingInd-1;
        if isempty(upperHomingInd)
            disp(['No homing achieved, subject ' subjectNames{f}(1) ' goodBeh ind ' num2str(i)])
            continue
        else
            homingTime(i) = interp1(distsFromEnd(lowerHomingInd:upperHomingInd),...
                kinTime(lowerHomingInd:upperHomingInd),homingEndDistance)...
                -midreachTime(i)-reactionTime(i)-goCueTime;
        end
        
        % If homing time is beyond trial limits + reaction
        % (homingStopAddedTime), exclude
        if homingTime(i)+midreachTime(i)+reactionTime(i) > maxHomingEndTime
            disp(['Homing too late, subject ' subjectNames{f}(1) ' goodBeh ind ' num2str(i), ', HaT+RT+HoT =  ' num2str(homingTime(i)+midreachTime(i)+reactionTime(i))])
            homingTime(i) = nan;
            continue
        end
    end; clear i
    
    % Arrange these by reward
    goodRewardLabels = [goodBeh.reward];
    for r = 1:nrewards
        curInds = goodRewardLabels==rewards(r);
        reachTime_bySubject_byReward{f,r} = homingTime(curInds)+midreachTime(curInds)+reactionTime(curInds);
        reactionTime_bySubject_byReward{f,r} = reactionTime(curInds);
        midreachTime_bySubject_byReward{f,r} = midreachTime(curInds);
        homingTime_bySubject_byReward{f,r} = homingTime(curInds);
    end; clear r
end; clear f

n_bySubject_byReward = cellfun(@(x) length(x), reactionTime_bySubject_byReward)
nBadHalfway_bySubject_byReward = cellfun(@(x) sum(isnan(x)), midreachTime_bySubject_byReward)
nBadHoming_bySubject_byReward = cellfun(@(x) sum(isnan(x)), homingTime_bySubject_byReward)

%% NOTE homing times are NOT normally distributed
qqplot(homingTime_bySubject_byReward{1,4})

%% as an example^ so mean isn't the best thing to use. Median is more
% appropriate. Let's get the medians and standard error of median for each
% group

nboots = 10000;
% homingTime_median_bySubject_byReward = cellfun(@(x) nanmedian(x), homingTime_bySubject_byReward);
homingTime_median_bySubject_byReward = nan(nsubjects,nrewards);
homingTime_semed_bySubject_byReward =  nan(nsubjects,nrewards);
for f = 1:nsubjects
    for r = 1:nrewards
            curVals = homingTime_bySubject_byReward{f,r};
            curVals(isnan(curVals))=[];
            n = length(curVals);
            bootMeds = nan(nboots,1);
            for b = 1:nboots
                bootMeds(b) = median(curVals(randsample(n,n,true))); % sample WITH replacement
            end; clear b
            homingTime_median_bySubject_byReward(f,r) = mean(bootMeds);
            homingTime_semed_bySubject_byReward(f,r) = std(bootMeds);
    end; clear r
end; clear f

%% Make a figure
rewNames = {'Small','Medium','Large','Jackpot'};
color = [0.5 0 0];
lw = 1;
ms = 1;
fs = 12;

figure
for f = 1:nsubjects
    subplot(1,nsubjects,f); hold on
    curMedian = homingTime_median_bySubject_byReward(f,:);
    curSemed = homingTime_semed_bySubject_byReward(f,:);
    errorbar(rewards, curMedian, curSemed,'k.-','linewidth',lw,'markersize',ms,'color',color);

    set(gca,'fontname','arial')
    set(gca,'fontsize',fs)
    xticks(1:nrewards)
    xticklabels(rewNames)
    xtickangle(45)
    yticks([90 150])
    if f == 1
        ylabel('Homing Time (ms)')
    else
        yticks([])
    end
    axis([rewards(1)-1 rewards(end)+1 87.5 152.5])
    title(['Monkey ' subjectNames{f}(1)])
end; clear f
set(gcf,'position',[33 638 895 319])

% Save it!
figname = 'Fig3C_medianHT';
saveas(gcf,['MATLABFigs\' figname])
saveas(gcf,['SVGs\' figname '.svg'])


%% Do stats - So many options, but I think the easiest thing to do is just
%  Mann-Whitney U-test for S->L and L->J. I know, I know, all the
%  statisticians will be like "bUt ThAt AkShuaLLy dOeSnT TeST tHe
%  DiffERenCe In MeDIaN!?!@11!1 *something something stochastically greater 
%  something something*" but for all intents and purposes it seems like it
%  works fine. If people are really upset, we'll just compare the
%  intervals found from bootstraps.
pvals_byRewardXReward_bySubject = nan(nrewards,nrewards,nsubjects);
for f = 1:nsubjects
    for r1 = 1:nrewards
        for r2 = 1:nrewards
            pvals_byRewardXReward_bySubject(r1,r2,f) = ranksum(homingTime_bySubject_byReward{f,r1},homingTime_bySubject_byReward{f,r2});
        end; clear r2
    end; clear r1
end; clear f

pvals_byRewardXReward_bySubject



%% display p values
for f = 1:nsubjects
    disp(['Subject ' subjectNames{f}(1)])
        disp(['S->L p = ' num2str(pvals_byRewardXReward_bySubject(1,3,f))])
        disp(['L->J p = ' num2str(pvals_byRewardXReward_bySubject(3,4,f))])
end; clear f

