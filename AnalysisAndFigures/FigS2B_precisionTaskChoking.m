%%% This script makes the figure to show the precision task performance.
%%%
%%% Adam Smoulder, 7/27/20 (edit 8/24/20)

% Load subject data behavior file
load('FordTubes_minimumTubeFeatures_27-Jul-2020133051');

%% Extract counts by reward, day, and status
trialStatusLabels = [behavior.trialStatusLabels];
statuses = unique(trialStatusLabels);
nstatuses = length(statuses);
rewardLabels = [behavior.reward];
rewards = unique(rewardLabels);
nrewards = length(rewards);
dayLabels = [behavior.day];
days = unique(dayLabels);
ndays = length(days);

n_byCond =  nan(nstatuses,nrewards,ndays);
for s = 1:nstatuses
    for r = 1:nrewards
        for a = 1:ndays
            n_byCond(s,r,a) = sum(trialStatusLabels==statuses(s) & ...
                rewardLabels==rewards(r) & ...
                dayLabels==days(a));
        end; clear a
    end; clear r
end; clear s

% Get success and total counts by reward
n_byReward = squeeze(sum(sum(n_byCond),3));
nsucc_byReward = squeeze(sum(n_byCond(2,:,:),3));


%% Do binomial proportion test and do bootstrapping
pvals_byRewardXReward = binomialProportionTest(nsucc_byReward,n_byReward);
nboots = 10000;
[~, nsucc_bootMean, ~, ~, ~, ~, nsucc_bootSem]...
    = bootstrapBinaryEvent(nsucc_byReward,n_byReward,nboots);
succRate_mean = nsucc_bootMean'./n_byReward*100;
succRate_sem = nsucc_bootSem'./n_byReward*100;

%% Get stuff for the last 8 sessions
nEndDaysToCheck = 8;
nend_byReward = squeeze(sum(sum(n_byCond(:,:,end-nEndDaysToCheck+1:end),3)));
nendsucc_byReward = squeeze(sum(n_byCond(2,:,end-nEndDaysToCheck+1:end),3));

pvalsend_byRewardXReward = binomialProportionTest(nendsucc_byReward,nend_byReward);
nboots = 10000;
[~, nendsucc_bootMean, ~, ~, ~, ~, nendsucc_bootSem]...
    = bootstrapBinaryEvent(nendsucc_byReward,nend_byReward,nboots);
succRateEnd_mean = nendsucc_bootMean'./nend_byReward*100;
succRateEnd_sem = nendsucc_bootSem'./nend_byReward*100;

%% Make a figure
rewNames = {'Small','Medium','Large','Jackpot'};

figure
errorbar(rewards-0.15,succRate_mean,succRate_sem,'k.-','linewidth',1,'markersize',20)
axis([0.5 0.5+nrewards 54 76])
set(gca,'fontname','arial')
set(gca,'fontsize',18)
xticks(1:nrewards)
xticklabels(rewNames)
% xtickangle(45)
yticks([55 65 75])
ylabel('Success Rate (%)')
title(['Monkey F, precision task ' num2str(ndays) ' sessions'])
set(gcf,'position',[418 548 521 400])

% on top of this, plot the last 8 sessions
hold on
errorbar(rewards+0.15,succRateEnd_mean,succRateEnd_sem,'k.-','linewidth',0.5,'markersize',1,'color',[128 128 0]/255)
% plot(rewards+0.1,succRateEnd_mean,'.-','linewidth',0.5,'markersize',20,'color',[0.5 0.5 0.5])
legend('All sessions', ['Last ' num2str(nEndDaysToCheck) ' sessions'],'location','SW')


pvals_byRewardXReward
n_byReward

pvalsend_byRewardXReward
nend_byReward

% Save it!
figname = 'FigS2B_precisionTaskChoking';
saveas(gcf,['MATLABFigs\' figname])
saveas(gcf,['SVGs\' figname '.svg'])
