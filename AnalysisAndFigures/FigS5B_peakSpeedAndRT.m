%%% This script makes the peak speed and reaction time panels for the
%%% idiosyncrasies figure
%%%
%%% Adam Smoulder, 8/18/20 (eidt 10/15/20 to add stats)

% Load subject data behavior files
loadBehaviorDataForMainFigures

% Get number of rewards; assume reward labels are same across subjects
rewards = unique([behavior_bySubject{1}.reward]);
nrewards = length(rewards);

%% Find the reaction (exit) time and peak speed by reward
homingStopDistance = 1; % mm - how far from end target to stop homing time calculation

reactionTime_bySubject_byReward = cell(nsubjects,nrewards);
peakSpeed_bySubject_byReward = cell(nsubjects,nrewards);

% Make a butterworth LPF for Earl's speed traces, since his haven't been
% filtered yet.
fsamp = 120;
fcut = 15;
order = 2; % using forward/backward filtering, so effectively 2*this = real order
[b,a] = butter(order,fcut/(fsamp/2));

for f = 1:nsubjects
    behavior = behavior_bySubject{f};
    startTargRad = behavior(1).task_parameters.startTargRadius;
    reachTargRad = behavior(1).task_parameters.endTargRadius;
    
    % Remove delay failures, quit outs, no-attempts
    goodBeh = behavior(~ismember([behavior.trialStatusLabels],[0 -111 -11 -12 -13 -14 -20 -21]));
    ngood = length(goodBeh);
    
    % Get all the timings
    reactionTime = nan(ngood,1);
    peakSpeed = nan(ngood,1);
    for i = 1:ngood
        % Get reaction time
        goCueTime = goodBeh(i).t3_goCueTimes;
        reactionTime(i) = goodBeh(i).t4c_reactionTimes_exit-goCueTime;
        
        % Get peak speed; we need to smooth the speed traces for Earl sadly
        if f ~=3
            peakSpeed(i) = goodBeh(i).peakSpeed;
        else
            curKin = goodBeh(i).kinematics_updated;
            smoothSpeed = filtfilt(b,a,curKin.speed);
            goodInds = curKin.time >= goCueTime & curKin.time <= goodBeh(i).t7_postReachStateTimes;
%             figure; hold on; plot(curKin.speed(goodInds)); plot(smoothSpeed(goodInds));
%             set(gcf,'position',[1109 49 560 420])
            peakSpeed(i) = max(smoothSpeed(goodInds));
        end
        
    end; clear i
    
    % Arrange these by reward
    goodRewardLabels = [goodBeh.reward];
    for r = 1:nrewards
        curInds = goodRewardLabels==rewards(r);
        reactionTime_bySubject_byReward{f,r} = reactionTime(curInds);
        peakSpeed_bySubject_byReward{f,r} = peakSpeed(curInds);
    end; clear r
end; clear f

n_bySubject_byReward = cellfun(@(x) length(x), reactionTime_bySubject_byReward)


%% Reaction times are not normal - hence we'll use the median. Since 
%  we're using the median, we need to find the empirical standard error of
%  the median with boostrapping...ugh
figure; qqplot(reactionTime_bySubject_byReward{1,1})

% as an example^ so mean isn't the best thing to use. Median is more
% appropriate. Let's get the medians and standard error of median for each
% group

nboots = 10000;
% homingTime_median_bySubject_byReward = cellfun(@(x) nanmedian(x), homingTime_bySubject_byReward);
reactionTime_median_bySubject_byReward = nan(nsubjects,nrewards);
reactionTime_semed_bySubject_byReward =  nan(nsubjects,nrewards);
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
    end; clear r
end; clear f

% Conversely, peak speed distrubtions are a bit more normal; I've tried 
% with median and the results hold anyways, so it doesn't really matter
figure; qqplot(peakSpeed_bySubject_byReward{1,1})
peakSpeed_mean_bySubject_byReward = cellfun(@(x) nanmean(x), peakSpeed_bySubject_byReward);
peakSpeed_median_bySubject_byReward = cellfun(@(x) nanmedian(x), peakSpeed_bySubject_byReward);
peakSpeed_sem_bySubject_byReward = cellfun(@(x) nansem(x,1), peakSpeed_bySubject_byReward);


%% Make a figure
rewNames = {'Small','Medium','Large','Jackpot'};
metricColors = [0 0 127; 255 0 255]/255;
lw = 1;
ms = 1;
fsamp = 10;

figure
for f = 1:nsubjects
    subplot(3,nsubjects,3+f); hold on
    curMedian = peakSpeed_mean_bySubject_byReward(f,:); 
    curSem = peakSpeed_sem_bySubject_byReward(f,:); 
    errorbar(rewards-0.2+0.1*f, curMedian, curSem,'-','linewidth',lw,'color',metricColors(1,:),'capsize',3);
    if f == 1
        ylabel('Peak Speed (m/s)')
    end
    yticks([0.39 0.48])
    set(gca,'fontsize',fsamp)
    set(gca,'fontname','arial')
    axis([rewards(1)-0.5 rewards(end)+0.5 .385 0.485])
    title(['Monkey ' subjectNames{f}(1)])
    
    
    subplot(3,nsubjects,6+f); hold on
    curMedian = reactionTime_median_bySubject_byReward(f,:); 
    curSem = reactionTime_semed_bySubject_byReward(f,:); 
    errorbar(rewards-0.2+0.1*f, curMedian, curSem,'-','linewidth',lw,'color',metricColors(2,:),'capsize',3);
    if f == 1
        ylabel('Reaction Time (ms)')
    end
    xticks([])
    yticks([225 290])
    set(gca,'fontsize',fsamp)
    axis([rewards(1)-0.5 rewards(end)+0.5 220 295])
    
    xticks(rewards);
    xticklabels(rewNames)
    xtickangle(45)
end; clear f
set(gcf,'position',[33 540 453 417])

% Save it!
figname = 'FigS5B_peakSpeedAndRT';
saveas(gcf,['MATLABFigs\' figname])
saveas(gcf,['SVGs\' figname '.svg'])



%% Do stats for each of them. First, peak speed. Since it's pretty normal,
%  we'll just use Welch's t-test
pvals_byRewardXReward_bySubject = nan(nrewards,nrewards,nsubjects);
for f = 1:nsubjects
    for r1 = 1:nrewards
        for r2 = 1:nrewards
            [~,pvals_byRewardXReward_bySubject(r1,r2,f),~,STATS] = ttest2(peakSpeed_bySubject_byReward{f,r1},peakSpeed_bySubject_byReward{f,r2});
        end; clear r2
    end; clear r1
end; clear f

pvals_byRewardXReward_bySubject

% display p values
for f = 1:nsubjects
    disp(['Subject ' subjectNames{f}(1)])
        disp(['S->L p = ' num2str(pvals_byRewardXReward_bySubject(1,3,f))])
        disp(['L->J p = ' num2str(pvals_byRewardXReward_bySubject(3,4,f))])
end; clear f


%% Now for reaction times. Since they're not super normal, we'll use a
%  Mann-Whitney U-test

pvals_byRewardXReward_bySubject = nan(nrewards,nrewards,nsubjects);
for f = 1:nsubjects
    for r1 = 1:nrewards
        for r2 = 1:nrewards
            pvals_byRewardXReward_bySubject(r1,r2,f) = ranksum(reactionTime_bySubject_byReward{f,r1},reactionTime_bySubject_byReward{f,r2});
        end; clear r2
    end; clear r1
end; clear f

pvals_byRewardXReward_bySubject

% display p values
for f = 1:nsubjects
    disp(['Subject ' subjectNames{f}(1)])
        disp(['S->L p = ' num2str(pvals_byRewardXReward_bySubject(1,3,f))])
        disp(['L->J p = ' num2str(pvals_byRewardXReward_bySubject(3,4,f))])
end; clear f





