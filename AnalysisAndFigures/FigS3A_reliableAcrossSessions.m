%%% This script makes the figure to show individual sessions' performance
%%% across rewards in the speed+accuracy task
%%%
%%% Adam Smoulder, 7/27/20

% Load subject data behavior files
loadBehaviorDataForMainFigures

%% Get success rates by day
excludeStatuses = [-21 -20 -13 -14 0]; % quitouts and misstarts are ignored in counts
n_bySubject_byDay_byReward = cell(nsubjects,1);
nsucc_bySubject_byDay_byReward = cell(nsubjects,1);
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
    
    % Populate the counts for each day x reward
    n_bySubject_byDay_byReward{f} = nan(ndays,nrewards);
    nsucc_bySubject_byDay_byReward{f} = nan(ndays,nrewards);
    for a = 1:ndays
        for r = 1:nrewards
            n_bySubject_byDay_byReward{f}(a,r) = ...
                sum(dayLabels==days(a) & rewardLabels==rewards(r) & ...
                ~ismember(trialStatusLabels,excludeStatuses));
            nsucc_bySubject_byDay_byReward{f}(a,r) = ...
                sum(dayLabels==days(a) & rewardLabels==rewards(r) & ...
                ismember(trialStatusLabels,1));   
        end; clear r
    end; clear a
end; clear f

% We calculate a success rate matrix for each subject
succRate_bySubject_byDay_byReward = cellfun(@(x,y) 100*x./y,...
    nsucc_bySubject_byDay_byReward, n_bySubject_byDay_byReward,...
    'UniformOutput',false);


%% Make the plot for each animal
rewNames = {'Small','Medium','Large','Jackpot'};
% colors = colormap(copper);
colors = ... % From a site of random distinguishable colors
[[230, 25, 75];
[60, 180, 75]; 
[255, 225, 25]; 
[0, 130, 200]; 
[245, 130, 48]; 
[70, 240, 240]; 
[240, 50, 230]; 
[0, 128, 128]; 
[220, 190, 255]; 
[170, 110, 40]; 
[128, 0, 0]; 
[170, 255, 195]; 
[0, 0, 128]; 
[128, 128, 128]; 
[255, 255, 255]; 
[0, 0, 0]]/255;

alpha = 1
lw = 0.5
ms = 10


figure
for f = 1:nsubjects
    subplot(1,nsubjects,f)
    hold on
    succRate = succRate_bySubject_byDay_byReward{f};
    ndays = size(succRate,1);
    curColors = colors(round(linspace(1,size(colors,1),ndays)),:);
    for a = 1:ndays
        plot(rewards+0.3/ndays*a-0.15,succRate(a,:),'k.-','color',[colors(a,:) alpha],'linewidth',lw,'markersize',ms)
    end; clear a
    axis([0.5 0.5+nrewards 17.5 92.5])
    set(gca,'fontname','arial')
    set(gca,'fontsize',16)
    xticks(1:nrewards)
    xticklabels(rewNames)
    xtickangle(45)
    yticks([20 55 90])
    if f == 1
        ylabel('Success Rate (%)')
    else
        yticks([])
    end
    title(['Monkey ' subjectNames{f}(1)])
    if f == 3
        legend('Day 1', 'Day 2', 'Day 3', 'Day 4','Day 5','Day 6','Day 7', 'Day 8','Day 9','Day 10','Day 11','location',[0.9 0.23 0.1 0.6])
    end
end; clear f
set(gcf,'position',[418 548 1042 400])


% Save it!
figname = 'FigS3A_reliableAcrossSessions';
saveas(gcf,['MATLABFigs\' figname])
saveas(gcf,['SVGs\' figname '.svg'])




