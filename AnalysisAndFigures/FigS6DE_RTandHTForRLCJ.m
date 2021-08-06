%%% This figure shows the reaction and homing time for the Rare-Large and
%%% Common-Jackpot controls.
%%%
%%% Adam Smoulder, 4/9/21

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


%% Get medians

nboots = 10000;
% homingTime_median_bySubject_byReward = cellfun(@(x) nanmedian(x), homingTime_bySubject_byReward);
reactionTime_median_bySubject_byReward = nan(nsubjects,nrewards);
reactionTime_semed_bySubject_byReward =  nan(nsubjects,nrewards);
homingTime_median_bySubject_byReward = nan(nsubjects,nrewards);
homingTime_semed_bySubject_byReward =  nan(nsubjects,nrewards);
for f = 1:nsubjects
    for r = 1:nrewards
            curVals = reactionTime_bySubject_byReward{f,r};
            curVals(isnan(curVals))=[];
            n = length(curVals);
            bootMeds = nan(nboots,1);
            for b = 1:nboots
                bootMeds(b) = median(curVals(randsample(n,n,true))); % sample WITH replacement
            end; clear b
            reactionTime_median_bySubject_byReward(f,r) = mean(bootMeds);
            reactionTime_semed_bySubject_byReward(f,r) = std(bootMeds);
        
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
RLColor = 0.75*[0 1 1];
CJColor = 0.9*[1 0 1];
RTColor = [0 0 0];
HTColor = [0 0 0];
lw = 1;
ms = 1;
fs = 12;

figure
for f = 1:nsubjects
    if f == 1, curColor = RLColor; else curColor = CJColor; end
    
    subplot(2,nsubjects,f); hold on
    curMedian = reactionTime_median_bySubject_byReward(f,1:4);
    curSemed = reactionTime_semed_bySubject_byReward(f,1:4);
    errorbar([1 2 3 3.92], curMedian, curSemed,'k.-','linewidth',lw,'markersize',ms,'color',curColor);
    curMedian = reactionTime_median_bySubject_byReward(f,[1 2 3 5]);
    curSemed = reactionTime_semed_bySubject_byReward(f,[1 2 3 5]);
    errorbar([1 2 3 4.08], curMedian, curSemed,'k.-','linewidth',lw,'markersize',ms,'color',RTColor);
    
    set(gca,'fontname','arial')
    set(gca,'fontsize',fs)
    set(gca,'TickDir','out')
    xticks([])
    
    xtickangle(45)
    yticks([230 310])
    if f == 1
        ylabel('Reaction Time (ms)')
        title('Rare-Large Sessions')
    else
        title('Common-Jackpot Sessions')
        yticks([])
    end
    axis([rewards(1)-0.5 4.5 225 315])
    
    
    subplot(2,nsubjects,f+nsubjects); hold on
    curMedian = homingTime_median_bySubject_byReward(f,1:4);
    curSemed = homingTime_semed_bySubject_byReward(f,1:4);
    errorbar([1 2 3 3.92], curMedian, curSemed,'k.-','linewidth',lw,'markersize',ms,'color',curColor);
    curMedian = homingTime_median_bySubject_byReward(f,[1 2 3 5]);
    curSemed = homingTime_semed_bySubject_byReward(f,[1 2 3 5]);
    errorbar([1 2 3 4.08], curMedian, curSemed,'k.-','linewidth',lw,'markersize',ms,'color',HTColor);
    
    set(gca,'fontname','arial')
    set(gca,'fontsize',fs)
    set(gca,'TickDir','out')
    xticks(1:4)
    if f == 1
        xticklabels({'S','M','L','RL/J'})
    else
        xticklabels({'S','M','L','CJ/J'})
    end
    
    xtickangle(45)
    yticks([90 135])
    if f == 1
        ylabel('Homing Time (ms)')
    else
        yticks([])
    end
    axis([rewards(1)-0.5 4.5 85 140])
end; clear f
set(gcf,'position',[33 458 648 499])

% Save it!
figname = 'FigS6DE_RTandHTForRLCJ';
saveas(gcf,['MATLABFigs\' figname])
saveas(gcf,['SVGs\' figname '.svg'])


%% Get p-vals for RT
pvalsRT_byRewardXReward_bySubject = nan(nrewards,nrewards,nsubjects);
for f = 1:nsubjects
    for r1 = 1:nrewards
        for r2 = 1:nrewards
            pvalsRT_byRewardXReward_bySubject(r1,r2,f) = ranksum(reactionTime_bySubject_byReward{f,r1},reactionTime_bySubject_byReward{f,r2});
        end; clear r2
    end; clear r1
end; clear f

pvalsRT_byRewardXReward_bySubject


%% Get p-vals for HT
pvalsHT_byRewardXReward_bySubject = nan(nrewards,nrewards,nsubjects);
for f = 1:nsubjects
    for r1 = 1:nrewards
        for r2 = 1:nrewards
            pvalsHT_byRewardXReward_bySubject(r1,r2,f) = ranksum(homingTime_bySubject_byReward{f,r1},homingTime_bySubject_byReward{f,r2});
        end; clear r2
    end; clear r1
end; clear f

pvalsHT_byRewardXReward_bySubject