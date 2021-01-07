%%% This script shows the breakdown of reach epoch statuses into
%%% overshoots, undershoots, and reach epoch success for each animal.
%%%
%%% Written by Patrick Marino, 10/18/20, edited by Adam Smoulder, 1/7/21

% Load subject data behavior files
loadBehaviorDataForMainFigures

% Get trials by direction, reward, and trial status conditions
getNCondFromBehaviorForMainFigures

% Only keep trials that make it to the reach; for reach epoch, we call
% anything that makes it to target hold (w/o overshooting) success (1), and 
% have overshoots (-22) and undershoots (-23)
for f = 1:nsubjects
    curBeh = behavior_bySubject{f};
    trialStatusLabels = trialStatusLabels_bySubject{f};
    trialStatusLabels(ismember(trialStatusLabels,[-22 -31 -32])) = -22;
    trialStatusLabels(ismember(trialStatusLabels,[-23 -24 -25])) = -23;
    trialsToKeep = ismember(trialStatusLabels,[1 -22 -23]);
    ngood = sum(trialsToKeep);
    goodBeh = curBeh(trialsToKeep);
    goodLabels = trialStatusLabels(trialsToKeep);
    
    % keep all the trajectories
    for i = 1:ngood
        goodBeh(i).trialStatusLabels = goodLabels(i);
    end; clear i
    behavior_bySubject{f} = goodBeh;
end; clear f

% Establish reward and status variables
rewards = 1:4;
rewardNames = {'Small','Medium','Large','Jackpot'};
numSubjects = length(subjectNames);
numRewards = length(rewardNames);
statusCodes = [1,-23,-22]; %Suc, Under, Over
numStatus = length(statusCodes);
%% Quick visual check to confirm each monkey's data only includes desired status codes
for subject = 1:numSubjects
    subjectData = behavior_bySubject{subject,1};
    subject
    uniqueTrialStatusLabels = unique([subjectData.trialStatusLabels])
end

%% Count up trials and trial statuses by subject and reward. Get status rates (as percentages)
totalReaches_bySubject_byReward = zeros(3,4);
statCounts_bySubject_byStatus_byReward = zeros(3,3,4);
statRates_bySubject_byStatus_byReward = zeros(3,3,4);
for subject = 1:numSubjects
    subjectData = behavior_bySubject{subject,1};
    for reward = 1:numRewards
        subjectRewardData = subjectData([subjectData.reward]==reward);
        totalReaches = size(subjectRewardData,2);
        totalReaches_bySubject_byReward(subject,reward) = totalReaches;
        for status = 1:numStatus
            statusCode = statusCodes(status);
            subjectRewardStatusData = subjectRewardData([subjectRewardData.trialStatusLabels]==statusCode);
            subjectRewardStatusCount = size(subjectRewardStatusData,2);
            statCounts_bySubject_byStatus_byReward(subject,status,reward) = subjectRewardStatusCount;
            statRates_bySubject_byStatus_byReward(subject,status,reward) = (subjectRewardStatusCount/totalReaches).*100;
        end
    end
end

%% Run binomial proportion test on status counts; Determine number of significance stars to include on final plot
bpStruct = struct('Subject',[],'Status',[],'pValues',[],'zStats',[],'hValues',[],'S_L_pVal',[],'S_L_numStars',[],'L_J_pVal',[],'L_J_numStars',[]);
structInd = 1;
for subject = 1:numSubjects
    %Get total number of reaches for each reward
    totalCount = totalReaches_bySubject_byReward(subject,:);
    for status = 1:numStatus
        %Run test for each status and store results
        eventCount = squeeze(statCounts_bySubject_byStatus_byReward(subject,status,:))';
        [pValues,zStats,hValues] = binomialProportionTest(eventCount,totalCount);
        bpStruct(structInd).Subject = subject;
        bpStruct(structInd).Status = statusCodes(status);
        bpStruct(structInd).pValues = pValues;
        bpStruct(structInd).zStats = zStats;
        bpStruct(structInd).hValues = hValues;
        %Get stars for S->L
        S_L_pVal = pValues(1,3);
        if S_L_pVal < 0.001
            numStars = 3;
        elseif S_L_pVal < 0.01
            numStars = 2;
        elseif S_L_pVal < 0.05
            numStars = 1;
        else
            numStars = 0;
        end
        bpStruct(structInd).S_L_pVal = S_L_pVal;
        bpStruct(structInd).S_L_numStars = numStars;
        %Get stars for L->J
        L_J_pVal = pValues(3,4);
        if L_J_pVal < 0.001
            numStars = 3;
        elseif L_J_pVal < 0.01
            numStars = 2;
        elseif L_J_pVal < 0.05
            numStars = 1;
        else
            numStars = 0;
        end
        bpStruct(structInd).L_J_pVal = L_J_pVal;
        bpStruct(structInd).L_J_numStars = numStars;
        structInd = structInd + 1;
    end
end

%% Create stacked bar plot of trial status rates
statusColors = [0 1 0; 1 1 1 ; 0 0 0 ] + [0 ; -0.2 ; 0.2];
lw = 1.5; %Line width
ms = 15; %Marker width
fs = 14; %Font size

figure('Position',[100 100 1250 400])
for subject = 1:numSubjects
    subplot(1,numSubjects,subject)
    %Get Fraction of Reaches (FOR) by status
    curSuccessFOR = squeeze(statRates_bySubject_byStatus_byReward(subject,1,:)); 
    curUshootFOR = squeeze(statRates_bySubject_byStatus_byReward(subject,2,:));
    curOshootFOR = squeeze(statRates_bySubject_byStatus_byReward(subject,3,:));
    %Bar Plot
    p = bar([curSuccessFOR curUshootFOR curOshootFOR],'stacked');
    p(1).FaceColor = statusColors(1,:);
    p(2).FaceColor = statusColors(2,:);
    p(3).FaceColor = statusColors(3,:);
    %Add Line
    hold on
    plot(rewards,curSuccessFOR,'k.-','linewidth',lw,'markersize',ms)
    %Format; add labels
    set(gca,'fontname','arial')
    set(gca,'fontsize',fs)
    axis([rewards(1)-1 rewards(end)+1 40 100])
    xticks(1:numRewards)
    xticklabels(rewardNames)
    xtickangle(45)
    yticks([40 50 100])
    yticklabels({'0', '50', '100'})
    if subject == 1
        ylabel('% of Reaches')
    end
     set(gca,'TickDir','out')
    title(['Monkey ' subjectNames{1,subject}(1)])
end; clear subject

% Save it!
figname = 'Fig4B_reachStatusBars';
saveas(gcf,['MATLABFigs\' figname])
saveas(gcf,['SVGs\' figname '.svg'])
