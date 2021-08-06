%%% This script makes the data figures showing ballistic reach endpoint
%%% predictions. 
%%%
%%% Adam Smoulder, 8/7/20

% Load subject data behavior files
loadBehaviorDataForMainFigures

% Get number of rewards; assume reward labels are same across subjects
rewards = unique([behavior_bySubject{1}.reward]);
nrewards = length(rewards);

BREP_byRew = cell(nsubjects,nrewards);
BREP_byRew_median = nan(nsubjects,nrewards,2);
BREP_byRew_mean = nan(nsubjects,nrewards,2);
BREP_byRew_semed = nan(nsubjects,nrewards,2);
BREP_byRew_sem = nan(nsubjects,nrewards,2);

%% Now really Looking at ballistic reach endpoint predictions

nboots = 10000; % for calculating SE of median

for f = 1:nsubjects
    % load the given name
    name = ['Monkey ' subjectNames{f}(1)];
    behavior = behavior_bySubject{f};
    
    % First, get out successes vs. failures
    trialStatusLabels = [behavior.trialStatusLabels];
   
    behaviorAll = behavior(~ismember([behavior.trialStatusLabels],[0 -111 -11 -12 -13 -14 -20 -21 -26]) & ...%     behaviorAll = behavior(ismember(trialStatusLabels,[1])' & ...
        ([behavior.t7_postReachStateTimes] > [behavior.t5_peakSpeedTimes]) & ...
        ([behavior.t5_peakSpeedTimes] > [behavior.t4_reactionTimes]));

    % For each reward, find the 3 values
    rewardLabels_all = [behaviorAll.reward];
    for r = 1:nrewards
        % For all data
        curInds = rewardLabels_all==r;
        curBeh = behaviorAll(curInds);
        curRotBREPs = nan(length(curBeh),2);
        for i = 1:length(curBeh)
            curRotBREPs(i,:) = curBeh(i).kinematics_updated.rotatedBallisticPredEndpoint;
        end; clear i
        BREP_byRew(f,r) = {curRotBREPs};
%         BREP_byRew_median(f,r,:) = median(curRotBREPs);
        BREP_byRew_mean(f,r,:) = mean(curRotBREPs);
        BREP_byRew_sem(f,r,:) = std(curRotBREPs)/sqrt(length(curRotBREPs));
        curMeds = nan(nboots,2);
        for b = 1:nboots
            curMeds(b,:) = median(curRotBREPs(randsample(length(curBeh),length(curBeh),true),:)); % sample WITH replacement
        end; clear b
        BREP_byRew_median(f,r,:) = mean(curMeds);
        BREP_byRew_semed(f,r,:) = std(curMeds);
    end; clear r
    disp(['Cleared ' num2str(f)' ' subjects of ' num2str(nsubjects)])
    
    endTargRadii(f) = behavior(1).task_parameters.endTargRadius;
end; clear n

%% Does the data seem normal?

figure; qqplot(curRotBREPs(:,1))
% yeah I'd say so.

%% and plot it!
figure; hold on
fs = 12;
% subplot(4,2,[7 8]); hold on

rewardNames = {'S','M','L','J'};
rewColors = {[1,0,0],[1,0.647,0],[0,0,1],[0,0,0]};
theta = linspace(0,2*pi,100);
for f = 1:nsubjects
    for r = 1:4
        errorbar(squeeze(BREP_byRew_mean(f,r,1)),squeeze(BREP_byRew_mean(f,r,2)),...
            squeeze(BREP_byRew_sem(f,r,2)),squeeze(BREP_byRew_sem(f,r,2)),...
            squeeze(BREP_byRew_sem(f,r,1)),squeeze(BREP_byRew_sem(f,r,1)),...
            '-','linewidth',1,'color',rewColors{r},'capsize',0);
    end; clear r
    axis([59 77 -3.75 3.75])
    xticks([60 75])

    yticks([-3.5 0 3.5])
    
    xlabel('On-Target Axis (mm)')
    ylabel('Off-Target Axis (mm)')
    set(gca,'fontname','arial')
    set(gca,'fontsize',fs)
end; clear f

set(gcf,'position',[39 617 650 319])

figname = 'Fig3B_BREPs';
saveas(gcf,['MATLABFigs\' figname])
saveas(gcf,['SVGs\' figname '.svg'])

%% Do stats. I guess we'll use Welch's t-test?
pvals_byRewardXReward_bySubject = nan(nrewards,nrewards,nsubjects);
for f = 1:nsubjects
    for r1 = 1:nrewards
        for r2 = 1:nrewards
            [~,pvals_byRewardXReward_bySubject(r1,r2,f),~,STATS] = ttest2(BREP_byRew{f,r1}(:,1),BREP_byRew{f,r2}(:,1));
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
