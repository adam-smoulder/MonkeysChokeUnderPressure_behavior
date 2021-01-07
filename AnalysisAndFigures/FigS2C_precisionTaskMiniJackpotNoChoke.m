%%% This script makes the figure showing that when a "mini-jackpot" reward
%%% size of 800uL was used in Monkey F's precision task, choking was very
%%% weak (insignificant)
%%%
%%% Adam Smoulder, 7/27/20

% Load subject data behavior file
load('DataFolder\FordTubesMJ_minimumTubeFeatures_27-Jul-2020133443');

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

%% Make a figure
rewNames = {'S','M','L','MJ'};

figure
errorbar(rewards,succRate_mean,succRate_sem,'k.-','linewidth',0.5,'markersize',20)
axis([0.5 0.5+nrewards 62.5 87.5])
set(gca,'fontname','arial')
set(gca,'fontsize',18)
xticks(1:nrewards)
xticklabels(rewNames)
% xtickangle(45)
yticks([65 75 85])
ylabel('Success Rate (%)')
title(['Monkey F, precision task, ' num2str(ndays) ' sessions'])
set(gcf,'position',[418 548 521 400])

pvals_byRewardXReward
n_byReward

% Save it!
figname = 'FigS2C_precisionTaskMiniJackpotNoChoke';
saveas(gcf,['MATLABFigs\' figname])
saveas(gcf,['SVGs\' figname '.svg'])

