%%% This script makes the figure showing performance for sessions where a
%%% "common-jackpot" cue was present; it showed up with the same frequency
%%% as S/M/L, but had the same magnitude as Jackpot, answering the question
%%% "can reward magnitude alone explain the performance drop that we see in 
%%% the main task?" - the answer is nah, though there might be a slight
%%% drop (though statistically insignificant).
%%%
%%% This one is actually very interesting for two reasons:
%%% 1) It means that in our setup, we need the rarity to provide enough
%%% pressure to induce choking. If he gets Jackpots willy-nilly, choking
%%% isn't consistent nor nearly as strong as with rare jackpots
%%% 2) He STILL choked for the normal (rare) Jackpots even with the
%%% Common-jackpots present. This is odd in the sense that they're worth
%%% equal amounts, just one happens to be rarer. So why is he choking? My
%%% bet is that there's some sort of conditioning that has occurred over
%%% his few months of seeing the normal Jackpot cue, such that now he's all
%%% psyched out about it. I'd bet you if you trained a new rare jackpot
%%% simultaneously with the common jackpot, he'd choke for neither. 
%%%
%%% Adam Smoulder, 7/27/20

% Load subject data behavior file
load('DataFolder\EarlCJ_behaviorFeatures_20210107');

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

% %% Make a figure
% rewNames = {'S','M','L','CJ','J'};
% 
% figure
% errorbar(rewards,succRate_mean,succRate_sem,'k.-','linewidth',0.5,'markersize',20)
% axis([0.5 0.5+nrewards 42.5 87.5])
% set(gca,'fontname','arial')
% set(gca,'fontsize',18)
% xticks(1:nrewards)
% xticklabels(rewNames)
% % xtickangle(45)
% yticks([45 65 85])
% ylabel('Success Rate (%)')
% title(['Monkey E, ' num2str(ndays) ' sessions'])
% set(gcf,'position',[418 548 521 400])
% 
% pvals_byRewardXReward
% n_byReward
% 
% 
% % Save it!
% figname = 'Fig3C_commonJackpotNoChoke';
% saveas(gcf,figname)
% saveas(gcf,[figname '.svg'])


%% Try making it a slightly diff style
CJColor = 0.9*[1 0 1];
rewNames = {'S','M','L','J (CJ)'};

figure; hold on
errorbar(rewards(1:4),succRate_mean(1:4),succRate_sem(1:4),'.-','linewidth',0.5,'markersize',20,'color',CJColor)
errorbar(rewards(1:4),succRate_mean([1 2 3 5]),succRate_sem([1 2 3 5]),'k.-','linewidth',0.5,'markersize',20)
axis([0.5 0.5+nrewards-1 42.5 87.5])
set(gca,'fontname','arial')
set(gca,'fontsize',18)
xticks(1:nrewards)
xticklabels(rewNames)
% xtickangle(45)
yticks([45 65 85])
ylabel('Success Rate (%)')
title(['Monkey E, ' num2str(ndays) ' sessions'])
set(gcf,'position',[418 548 521 400])

pvals_byRewardXReward
n_byReward

% Save it!
figname = 'Fig2B_commonJackpotNoChoke';
saveas(gcf,['MATLABFigs\' figname])
saveas(gcf,['SVGs\' figname '.svg'])

%% display p values and success rate differences
disp(['S->L, ' num2str(diff(succRate_mean([1 3]))) '% p = ' num2str(pvals_byRewardXReward(1,3))])
disp(['L->CJ, ' num2str(diff(succRate_mean([3 4]))) '% p = ' num2str(pvals_byRewardXReward(4,3))])
disp(['L->J, ' num2str(diff(succRate_mean([3 5]))) '% p = ' num2str(pvals_byRewardXReward(5,3))])
