%%% This script makes the figure showing performance for sessions where a
%%% "rare-large" cue was pressent; it showed up with the same frequency as
%%% Jackpots, but had the same magnitude as Large, answering the question
%%% "does rarity (and perhaps surprise) alone explain the performance drop
%%% that we see in the main task?" - the answer is nah, though there might
%%% be a slight drop (albeit insignificant)
%%%
%%% Adam Smoulder, 7/27/20

% Load subject data behavior file
load('EarlRL_behaviorFeatures_20210107');

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

% Get success and total counts by reward (skip bad statuses)
badStatuses = [0 -111 -13 -14 -20 -21 -33];
n_byReward = squeeze(sum(sum(n_byCond(~ismember(statuses,badStatuses),:,:),3)));
nsucc_byReward = squeeze(sum(n_byCond(find(statuses==1,1),:,:),3));


%% Do binomial proportion test and do bootstrapping
pvals_byRewardXReward = binomialProportionTest(nsucc_byReward,n_byReward);
nboots = 10000;
[~, nsucc_bootMean, ~, ~, ~, ~, nsucc_bootSem]...
    = bootstrapBinaryEvent(nsucc_byReward,n_byReward,nboots);
succRate_mean = nsucc_bootMean'./n_byReward*100;
succRate_sem = nsucc_bootSem'./n_byReward*100;


%% Make figure
RLColor = 0.75*[0 1 1];
rewNames = {'S','M','L','J (RL)'};

figure; hold on
errorbar(rewards(1:4),succRate_mean(1:4),succRate_sem(1:4),'.-','linewidth',0.5,'markersize',20,'color',RLColor)
errorbar(rewards(1:4),succRate_mean([1 2 3 5]),succRate_sem([1 2 3 5]),'k.-','linewidth',0.5,'markersize',20)
axis([0.5 0.5+nrewards-1 47.5 92.5])
set(gca,'fontname','arial')
set(gca,'fontsize',18)
xticks(1:nrewards)
xticklabels(rewNames)
% xtickangle(45)
yticks([50 70 90])
ylabel('Success Rate (%)')
title(['Monkey E, ' num2str(ndays) ' sessions'])
set(gcf,'position',[418 548 521 400])

pvals_byRewardXReward
n_byReward

% Save it!
figname = 'Fig2A_rareLargeNoChoke';
saveas(gcf,['MATLABFigs\' figname])
saveas(gcf,['SVGs\' figname '.svg'])

%% display p values and differences in success rates

disp(['S->L, ' num2str(diff(succRate_mean([1 3]))) '% p = ' num2str(pvals_byRewardXReward(1,3))])
disp(['L->RL, ' num2str(diff(succRate_mean([3 4]))) '% p = ' num2str(pvals_byRewardXReward(4,3))])
disp(['L->J, ' num2str(diff(succRate_mean([3 5]))) '% p = ' num2str(pvals_byRewardXReward(5,3))])

