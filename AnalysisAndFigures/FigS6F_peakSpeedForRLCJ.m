%%% This script makes the peak speed figure for the control experiments
%%%
%%% Adam Smoulder, 4/12/21

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

% Get number of rewards; assume reward labels are same across subjects
rewards = unique([behavior_bySubject{1}.reward]);
nrewards = length(rewards);

%% Find the reaction (exit) time and peak speed by reward
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

% Conversely, peak speed distrubtions are rather normal
peakSpeed_mean_bySubject_byReward = cellfun(@(x) nanmean(x), peakSpeed_bySubject_byReward);
peakSpeed_median_bySubject_byReward = cellfun(@(x) nanmedian(x), peakSpeed_bySubject_byReward);
peakSpeed_sem_bySubject_byReward = cellfun(@(x) nansem(x,1), peakSpeed_bySubject_byReward);

%% Make a figure
rewNames = {'Small','Medium','Large','Jackpot'};
RLColor = 0.75*[0 1 1];
CJColor = 0.9*[1 0 1];
metricColors = [0 0 0];
lw = 1;
ms = 1;
fsamp = 10;

figure
for f = 1:nsubjects
    if f == 1, curColor = RLColor; else curColor = CJColor; end
    
    subplot(1,nsubjects,f); hold on
    curMedian = peakSpeed_mean_bySubject_byReward(f,1:4); 
    curSem = peakSpeed_sem_bySubject_byReward(f,1:4); 
    errorbar([1 2 3 3.92], curMedian, curSem,'-','linewidth',lw,'color',curColor,'capsize',3);
    curMedian = peakSpeed_mean_bySubject_byReward(f,[1 2 3 5]); 
    curSem = peakSpeed_sem_bySubject_byReward(f,[1 2 3 5]); 
    errorbar([1 2 3 4.08], curMedian, curSem,'-','linewidth',lw,'color',metricColors(1,:),'capsize',3);
    
    if f == 1
        ylabel('Peak Speed (m/s)')
    end
    yticks([0.39 0.54])
    set(gca,'fontsize',fsamp)
    set(gca,'fontname','arial')
    axis([rewards(1)-0.5 4.5 .385 0.545])
    if f ==1
        title('Rare-Large Sessions')
    else
        title('Common-Jackpot Sessions')
    end
    
    if f == 1
        xticklabels({'S','M','L','RL/J'})
    else
        xticklabels({'S','M','L','CJ/J'})
    end
    
end; clear f
set(gcf,'position',[598 488 648 211])

% Save it!
figname = 'FigS6F_peakSpeedForRLCJ';
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






