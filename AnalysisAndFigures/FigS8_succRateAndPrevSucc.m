%%% This script looks at how previous trials failing affects success rates
%%%
%%% Adam Smoulder, 10/30/2020 (edit 5/4/21)

% Load subject data behavior files
loadBehaviorDataForMainFigures

% Get trials by direction, reward, and trial status conditions
getNCondFromBehaviorForMainFigures

%% Q1) Look at current trial; given previous trial succ/fail, outcome of
%  current trial?
n_bySubject_byReward_byPrevSF_bySF = nan([nsubjects,nrewards,2,2]);
for f = 1:nsubjects
    trialStatusLabels = trialStatusLabels_bySubject{f};
    rewardLabels = rewardLabels_bySubject{f};
    for r = 1:nrewards
        % Get current indices, previous trial status, and present status
        curInds = find(rewardLabels==rewards(r));
        curInds(curInds==1) = []; % can't use first index; no previous
        prevTrialSF = trialStatusLabels(curInds-1)==1;
        presTrialSF = trialStatusLabels(curInds)==1;
        n_bySubject_byReward_byPrevSF_bySF(f,r,1,1) = sum(prevTrialSF & presTrialSF);   % prev succ, cur succ
        n_bySubject_byReward_byPrevSF_bySF(f,r,1,2) = sum(prevTrialSF & ~presTrialSF);  % prev succ, cur fail
        n_bySubject_byReward_byPrevSF_bySF(f,r,2,1) = sum(~prevTrialSF & presTrialSF);  % prev fail, cur succ
        n_bySubject_byReward_byPrevSF_bySF(f,r,2,2) = sum(~prevTrialSF & ~presTrialSF); % prev fail, cur fail
    end; clear r
end; clear f

% Get rates for within each reward/pres condition
succRate_bySubject_byReward_byPresSF = ...
    100*n_bySubject_byReward_byPrevSF_bySF(:,:,:,1)./sum(n_bySubject_byReward_byPrevSF_bySF,4);

%% Do bootstrapping and binomial proportion test based on these counts
nboots = 10000;
succRate_mean_bySubject_byReward_byPrevSF = nan(nsubjects,nrewards,2);
succRate_sem_bySubject_byReward_byPrevSF = nan(nsubjects,nrewards,2);
pvals_byRewardXReward_bySubject_byPrevSF = nan(nrewards,nrewards,nsubjects,2);
for f = 1:nsubjects
    for p = 1:2 % 1 = prev succ, 2 = prev fail
        eventCounts = squeeze(n_bySubject_byReward_byPrevSF_bySF(f,:,p,1)); % successes are the first page (dim 3); we marginalize over directions
        totalCounts = squeeze(sum(n_bySubject_byReward_byPrevSF_bySF(f,:,p,:),4)); % all trials by reward size
        [~, bootMean, ~, ~, bar95low, bar95up, bootSem] = ...
            bootstrapBinaryEvent(eventCounts,totalCounts,nboots);
        succRate_mean_bySubject_byReward_byPrevSF(f,:,p) = bootMean./totalCounts'*100;
        succRate_sem_bySubject_byReward_byPrevSF(f,:,p) = bootSem./totalCounts'*100;
        
        % Use the binomial proportion test for significance; not related to
        % bootstrap, but might as well nab it here
        pvals_byRewardXReward_bySubject_byPrevSF(:,:,f,p) = binomialProportionTest(eventCounts,totalCounts);
    end; clear p
    disp(['Completed bootstrapping success rates for ' subjectNames{f}])
end; clear f

%% Plot it!
rewNames = {'Small','Medium','Large','Jackpot'};
prevColors = {[0 0 0.8],[0.8 0 0]};

figure
for f = 1:nsubjects
    subplot(1,nsubjects,f)
    hold on
    for p = 2:-1:1 % 1 = prev succ, 2 = prev fail; we do backwards to order red then blue
        errorbar(rewards-0.1+(p-1)*0.2,squeeze(succRate_mean_bySubject_byReward_byPrevSF(f,:,p)),...
            squeeze(succRate_sem_bySubject_byReward_byPrevSF(f,:,p)),...
            'k.-','linewidth',0.5,'markersize',20,'color',prevColors{p})
    end; clear p
    if f == 1
        legend('Prev trial fail','Prev trial succ', 'location','sw')
    end
    axis([0.5 0.5+nrewards 30 90])
    set(gca,'fontname','arial')
    set(gca,'fontsize',16)
    xticks(1:nrewards)
    xticklabels(rewNames)
    xtickangle(45)
    yticks([35 60 85])
    if f == 1
        ylabel('Success Rate (%)')
    end
    if f == 2
        xlabel('Reward')
    end
    title(['Monkey ' subjectNames{f}(1) ', ' num2str(length(unique([behavior_bySubject{f}.day]))) ' sessions'])
end; clear f
set(gcf,'position',[418 548 1042 400])


pvals_byRewardXReward_bySubject_byPrevSF

% Save it!
figname = 'FigS8_succRateAndPrevSucc';
saveas(gcf,['MATLABFigs\' figname])
saveas(gcf,['SVGs\' figname '.svg'])




% %% Q2) the reverse; disregard the current reward label. Based on the
% %  previous reward label and it's status, what are the odds of success?
% 
% n_bySubject_byPrevReward_byPrevSF_bySF = nan([nsubjects,nrewards,2,2]);
% for f = 1:nsubjects
%     trialStatusLabels = trialStatusLabels_bySubject{f};
%     rewardLabels = rewardLabels_bySubject{f};
%     prevRewardLabels = [0 ; rewardLabels(1:end-1)];
%     for r = 1:nrewards
%         % Get current indices, previous trial status, and present status
%         curInds = find(prevRewardLabels==rewards(r));
%         curInds(curInds==1) = []; % can't use first index; no previous
%         prevTrialSF = trialStatusLabels(curInds-1)==1;
%         presTrialSF = trialStatusLabels(curInds)==1;
%         n_bySubject_byPrevReward_byPrevSF_bySF(f,r,1,1) = sum(prevTrialSF & presTrialSF);   % prev succ, cur succ
%         n_bySubject_byPrevReward_byPrevSF_bySF(f,r,1,2) = sum(prevTrialSF & ~presTrialSF);  % prev succ, cur fail
%         n_bySubject_byPrevReward_byPrevSF_bySF(f,r,2,1) = sum(~prevTrialSF & presTrialSF);  % prev fail, cur succ
%         n_bySubject_byPrevReward_byPrevSF_bySF(f,r,2,2) = sum(~prevTrialSF & ~presTrialSF); % prev fail, cur fail
%     end; clear r
% end; clear f
% 
% % Get rates for within each reward/pres condition
% succRate_bySubject_byPrevReward_byPresSF = ...
%     100*n_bySubject_byPrevReward_byPrevSF_bySF(:,:,:,1)./sum(n_bySubject_byPrevReward_byPrevSF_bySF,4);
% 
% %% Do bootstrapping and binomial proportion test based on these counts
% nboots = 10000;
% succRate_mean_bySubject_byReward_byPrevSF = nan(nsubjects,nrewards,2);
% succRate_sem_bySubject_byReward_byPrevSF = nan(nsubjects,nrewards,2);
% pvals_byPrevRewardXPrevReward_bySubject_byPrevSF = nan(nrewards,nrewards,nsubjects,2);
% for f = 1:nsubjects
%     for p = 1:2 % 1 = prev succ, 2 = prev fail
%         eventCounts = squeeze(n_bySubject_byPrevReward_byPrevSF_bySF(f,:,p,1)); % successes are the first page (dim 3); we marginalize over directions
%         totalCounts = squeeze(sum(n_bySubject_byPrevReward_byPrevSF_bySF(f,:,p,:),4)); % all trials by reward size
%         [~, bootMean, ~, ~, bar95low, bar95up, bootSem] = ...
%             bootstrapBinaryEvent(eventCounts,totalCounts,nboots);
%         succRate_mean_bySubject_byReward_byPrevSF(f,:,p) = bootMean./totalCounts'*100;
%         succRate_sem_bySubject_byReward_byPrevSF(f,:,p) = bootSem./totalCounts'*100;
%         
%         % Use the binomial proportion test for significance; not related to
%         % bootstrap, but might as well nab it here
%         pvals_byPrevRewardXPrevReward_bySubject_byPrevSF(:,:,f,p) = binomialProportionTest(eventCounts,totalCounts);
%     end; clear p
%     disp(['Completed bootstrapping success rates for ' subjectNames{f}])
% end; clear f
% 
% 
% %% Plot it!
% rewNames = {'Small','Medium','Large','Jackpot'};
% prevColors = {[0 0.8 0],[0.8 0 0]};
% 
% figure
% for f = 1:nsubjects
%     subplot(1,nsubjects,f)
%     hold on
%     for p = 1:2 % 1 = prev succ, 2 = prev fail
%         errorbar(rewards-0.1+(p-1)*0.2,squeeze(succRate_mean_bySubject_byReward_byPrevSF(f,:,p)),...
%             squeeze(succRate_sem_bySubject_byReward_byPrevSF(f,:,p)),...
%             'k.-','linewidth',0.5,'markersize',20,'color',prevColors{p})
%     end; clear p
%     if f == 1
%         legend('Prev trial succ','Prev trial fail', 'location','sw')
%     end
%     axis([0.5 0.5+nrewards 30 90])
%     set(gca,'fontname','arial')
%     set(gca,'fontsize',16)
%     xticks(1:nrewards)
%     xticklabels(rewNames)
%     xtickangle(45)
%     yticks([35 60 85])
%     if f == 1
%         ylabel('Success Rate (%)')
%     end
%     if f == 2
%         xlabel('Previous trial reward size')
%     end
%     title(['Monkey ' subjectNames{f}(1) ', ' num2str(length(unique([behavior_bySubject{f}.day]))) ' sessions'])
% end; clear f
% set(gcf,'position',[418 548 1042 400])
% 
% 
% pvals_byPrevRewardXPrevReward_bySubject_byPrevSF
% 
% % Save it!
% figname = 'FigX_succRateAndPrevSuccByPrevRew';
% saveas(gcf,[figname])
% saveas(gcf,[figname '.svg'])



% %% Q3) Just show the full 4x4 for each animal.
% n_bySubject_byReward_byPrevRew_byPrevSF_bySF = nan([nsubjects,nrewards,nrewards,2,2]);
% for f = 1:nsubjects
%     trialStatusLabels = trialStatusLabels_bySubject{f};
%     rewardLabels = rewardLabels_bySubject{f};
%     prevRewardLabels = [0; rewardLabels(1:end-1)];
%     for r = 1:nrewards
%         for pr = 1:nrewards
%             % Get current indices, previous trial status, and present status
%             curInds = find(rewardLabels==rewards(r) & prevRewardLabels==rewards(pr));
%             prevTrialSF = trialStatusLabels(curInds-1)==1;
%             presTrialSF = trialStatusLabels(curInds)==1;
%             n_bySubject_byReward_byPrevRew_byPrevSF_bySF(f,r,pr,1,1) = sum(prevTrialSF & presTrialSF);   % prev succ, cur succ
%             n_bySubject_byReward_byPrevRew_byPrevSF_bySF(f,r,pr,1,2) = sum(prevTrialSF & ~presTrialSF);  % prev succ, cur fail
%             n_bySubject_byReward_byPrevRew_byPrevSF_bySF(f,r,pr,2,1) = sum(~prevTrialSF & presTrialSF);  % prev fail, cur succ
%             n_bySubject_byReward_byPrevRew_byPrevSF_bySF(f,r,pr,2,2) = sum(~prevTrialSF & ~presTrialSF); % prev fail, cur fail
%         end; clear pr
%     end; clear r
% end; clear f
% 
% % Get rates for within each reward/pres condition
% succRate_bySubject_byReward_byPrevRew_byPrevSF = ...
%     100*n_bySubject_byReward_byPrevRew_byPrevSF_bySF(:,:,:,:,1)./sum(n_bySubject_byReward_byPrevRew_byPrevSF_bySF,5);
% 
% 
% %% Stats and error bars
% nboots = 10000;
% succRate_mean_bySubject_byReward_byPrevRew_byPrevSF = nan(nsubjects,nrewards,nrewards,2);
% succRate_sem_bySubject_byReward_byPrevRew_byPrevSF = nan(nsubjects,nrewards,nrewards,2);
% pvals_byRewardXReward_bySubject_byPrevRew_byPrevSF = nan(nrewards,nrewards,nrewards,nsubjects,2);
% for f = 1:nsubjects
%     for p = 1:2 % 1 = prev succ, 2 = prev fail
%         for pr = 1:nrewards
%             eventCounts = squeeze(n_bySubject_byReward_byPrevRew_byPrevSF_bySF(f,:,pr,p,1)); % successes are the first page (dim 3); we marginalize over directions
%             totalCounts = squeeze(sum(n_bySubject_byReward_byPrevRew_byPrevSF_bySF(f,:,pr,p,:),5)); % all trials by reward size
%             [~, bootMean, ~, ~, bar95low, bar95up, bootSem] = ...
%                 bootstrapBinaryEvent(eventCounts,totalCounts,nboots);
%             succRate_mean_bySubject_byReward_byPrevRew_byPrevSF(f,:,pr,p) = bootMean./totalCounts'*100;
%             succRate_sem_bySubject_byReward_byPrevRew_byPrevSF(f,:,pr,p) = bootSem./totalCounts'*100;
%             
%             % Use the binomial proportion test for significance; not related to
%             % bootstrap, but might as well nab it here
%             pvals_byRewardXReward_bySubject_byPrevRew_byPrevSF(:,:,f,pr,p) = binomialProportionTest(eventCounts,totalCounts);
%         end; clear pr
%     end; clear p
%     disp(['Completed bootstrapping success rates for ' subjectNames{f}])
% end; clear f
% 
% 
% %% Plot it
% rewNames = {'Small','Medium','Large','Jackpot'};
% rewColors = {[1,0,0],[1,0.647,0],[0,0,1],[0,0,0]};
% figure;
% for f = 1:nsubjects
%     for p = 1:2
%         curData = squeeze(succRate_bySubject_byReward_byPrevRew_byPrevSF(f,:,:,p));
%         subplot(2,3,f+(p-1)*3); hold on
%         for pr = 1:nrewards
%             %             plot(curData(:,pr),'.-','linewidth',2,'markersize',20,'color',rewColors{pr});
%             errorbar(rewards-0.35+(pr)*0.1,squeeze(succRate_mean_bySubject_byReward_byPrevRew_byPrevSF(f,:,pr,p)),...
%                 squeeze(succRate_sem_bySubject_byReward_byPrevRew_byPrevSF(f,:,pr,p)),...
%                 'k.-','linewidth',0.5,'markersize',20,'color',rewColors{pr})
%         end; clear pr
%         axis([0.5 nrewards+0.5 0 100])
%         yticks([0 25 50 75 100])
%         set(gca,'fontname','arial')
%         set(gca,'fontsize',16)
%         xticks(1:nrewards)
%         xticklabels(rewNames)
%         xtickangle(30)
%         if f == 1 && p == 1
%             ylabel('Success Rate (%), prev trial SUCC')
%             legend('Prev S','Prev M','Prev L','Prev J','location','SE')
%         end
%         if f == 1 && p == 2
%             ylabel('Success Rate (%), prev trial FAIL')
%         end
%         title(['Monkey ' subjectNames{f}(1) ', ' num2str(length(unique([behavior_bySubject{f}.day]))) ' sessions'])
%     end; clear p
% end; clear f
% 
% set(gcf,'position',[1 41 1680 933])
% 
% % Save it!
% figname = 'FigX_succRateAndPreviousTrial_full';
% saveas(gcf,[figname])
% saveas(gcf,[figname '.svg'])
% 
% 
% 
% %% Ignoring SF of previous trial
% n_bySubject_byReward_byPrevRew_bySF = squeeze(sum(n_bySubject_byReward_byPrevRew_byPrevSF_bySF,4));
% 
% succRate_bySubject_byReward_byPrevRew = ...
%     100*n_bySubject_byReward_byPrevRew_bySF(:,:,:,1)./sum(n_bySubject_byReward_byPrevRew_bySF,4);
% 
% nboots = 10000;
% succRate_mean_bySubject_byReward_byPrevRew = nan(nsubjects,nrewards,nrewards);
% succRate_sem_bySubject_byReward_byPrevRew = nan(nsubjects,nrewards,nrewards);
% pvals_byRewardXReward_bySubject_byPrevRew = nan(nrewards,nrewards,nrewards,nsubjects);
% for f = 1:nsubjects
%     for pr = 1:nrewards
%         eventCounts = squeeze(n_bySubject_byReward_byPrevRew_bySF(f,:,pr,1)); % successes are the first page (dim 3); we marginalize over directions
%         totalCounts = squeeze(sum(n_bySubject_byReward_byPrevRew_bySF(f,:,pr,:),4)); % all trials by reward size
%         [~, bootMean, ~, ~, bar95low, bar95up, bootSem] = ...
%             bootstrapBinaryEvent(eventCounts,totalCounts,nboots);
%         succRate_mean_bySubject_byReward_byPrevRew(f,:,pr) = bootMean./totalCounts'*100;
%         succRate_sem_bySubject_byReward_byPrevRew(f,:,pr) = bootSem./totalCounts'*100;
%         
%         % Use the binomial proportion test for significance; not related to
%         % bootstrap, but might as well nab it here
%         pvals_byRewardXReward_bySubject_byPrevRew(:,:,f,pr) = binomialProportionTest(eventCounts,totalCounts);
%     end; clear pr
%     disp(['Completed bootstrapping success rates for ' subjectNames{f}])
% end; clear f
% 
% 
% %% Plot it
% 
% rewNames = {'Small','Medium','Large','Jackpot'};
% rewColors = getDistinctColors('SELECT_ORDER',8);
% figure;
% for f = 1:nsubjects
%     curData = squeeze(succRate_bySubject_byReward_byPrevRew(f,:,:));
%     subplot(1,3,f); hold on
%     for pr = 1:nrewards
%         errorbar(rewards-0.35+(pr)*0.1,squeeze(succRate_mean_bySubject_byReward_byPrevRew(f,:,pr)),...
%             squeeze(succRate_sem_bySubject_byReward_byPrevRew(f,:,pr)),...
%             'k.-','linewidth',0.5,'markersize',20,'color',rewColors{pr})
%     end; clear pr
%     axis([0.5 nrewards+0.5 15 90])
%     yticks([20 40 60 80])
%     set(gca,'fontname','arial')
%     set(gca,'fontsize',16)
%     xticks(1:nrewards)
%     xticklabels(rewNames)
%     xtickangle(45)
%     if f == 1
%         ylabel('Success Rate (%)')
%         legend('Prev S','Prev M','Prev L','Prev J','location','SW')
%     end
%     if f == 1
%         ylabel('Success Rate (%)')
%     end
%     if f == 2
%         xlabel('Reward')
%     end
%     title(['Monkey ' subjectNames{f}(1) ', ' num2str(length(unique([behavior_bySubject{f}.day]))) ' sessions'])
% end; clear f
% 
% set(gcf,'position',[-1828 409 1066 443])
% 
% % Save it!
% figname = 'FigX_succRateByPrevRew';
% saveas(gcf,[figname])
% saveas(gcf,[figname '.svg'])
