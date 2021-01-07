%%% This script makes the figure to show the average performance during the
%%% first half of sessions and the second half of sessions for each
%%% subject.
%%%
%%% Adam Smoulder, 7/27/20

% Load subject data behavior files
loadBehaviorDataForMainFigures

%% Get success rates by day
excludeStatuses = [-20 -13 -14 0]; % quitouts and misstarts are ignored in counts
n_bySubject_byDay_byHalf_byReward = cell(nsubjects,1);
nsucc_bySubject_byDay_byHalf_byReward = cell(nsubjects,1);
for f = 1:nsubjects
    % Get the current behavior, # of days, etc.
    behavior = behavior_bySubject{f};
    dayLabels = [behavior.day];
    days = unique(dayLabels);
    ndays = length(days);
    rewardLabels = [behavior.reward];
    rewards = unique(rewardLabels);
    nrewards = length(rewards);
    trialStatusLabels = [behavior.trialStatusLabels];
    
    % Make "half labels" for each day based on whether each trial is within 
    % first or second half of the session.
    halfLabels = nan(1,length(trialStatusLabels));
    for a = 1:ndays
        curTrials = find(dayLabels==a);
        halfLabels(curTrials(1:floor(length(curTrials)/2))) = 1;
        halfLabels(curTrials(floor(length(curTrials)/2)+1:end)) = 2;
    end; clear a
    
    
    % Populate the counts for each day x half xreward
    n_bySubject_byDay_byHalf_byReward{f} = nan(ndays,2,nrewards);
    nsucc_bySubject_byDay_byHalf_byReward{f} = nan(ndays,2,nrewards);
    for a = 1:ndays
        for r = 1:nrewards
            for h = 1:2
                n_bySubject_byDay_byHalf_byReward{f}(a,h,r) = ...
                    sum(dayLabels==days(a) & rewardLabels==rewards(r) & ...
                    ~ismember(trialStatusLabels,excludeStatuses) & ...
                    halfLabels==h);
                nsucc_bySubject_byDay_byHalf_byReward{f}(a,h,r) = ...
                    sum(dayLabels==days(a) & rewardLabels==rewards(r) & ...
                    ismember(trialStatusLabels,1) & ...
                    halfLabels==h);
            end; clear h
        end; clear r
    end; clear a
end; clear f

% We calculate a success rate matrix for each subject
succRate_bySubject_byHalf_byReward = cellfun(@(x,y) 100*squeeze(sum(x)./sum(y)),...
    nsucc_bySubject_byDay_byHalf_byReward, n_bySubject_byDay_byHalf_byReward,...
    'UniformOutput',false);


%% Calculate significance via binomial proportions test
pvals_byRewardXReward_bySubject_byHalf = nan(nrewards,nrewards,nsubjects,2);
for f = 1:nsubjects
    for h = 1:2
        pvals_byRewardXReward_bySubject_byHalf(:,:,f,h) = binomialProportionTest(...
            squeeze(sum(nsucc_bySubject_byDay_byHalf_byReward{f}(:,h,:)))',...
            squeeze(sum(n_bySubject_byDay_byHalf_byReward{f}(:,h,:)))');
    end; clear h
end; clear f


%% Make the plot for each animal
rewNames = {'Small','Medium','Large','Jackpot'};
subjectColors = [0 100 0 ; 0 180 0 ; 150 200 0]/255;
alpha = 1
lw = 2
ms = 10

figure
for f = 1:nsubjects
    subplot(1,nsubjects,f)
    hold on
    succRates = succRate_bySubject_byHalf_byReward{f};
    curColor = subjectColors(f,:);
    plot(rewards,succRates(1,:),'k^-','color',[curColor alpha],'linewidth',lw,'markersize',ms) % 1st half
    plot(rewards,succRates(2,:),'ks--','color',[curColor alpha],'linewidth',lw,'markersize',ms) % 2nd half
    axis([0.5 0.5+nrewards 37.5 92.5])
    set(gca,'fontname','arial')
    set(gca,'fontsize',16)
    xticks(1:nrewards)
    xticklabels(rewNames)
    xtickangle(45)
    yticks([40 65 90])
    if f == 1
        ylabel('Success Rate (%)')
    end
    if f == 3
        legend('1^{st} Half','2^{nd} Half','location','SE')
    end
    title(['Monkey ' subjectNames{f}(1)])
end; clear f
set(gcf,'position',[418 548 1042 400])


% Save it!
figname = 'FigS3B_reliableWithinSessions';
saveas(gcf,['MATLABFigs\' figname])
saveas(gcf,['SVGs\' figname '.svg'])
