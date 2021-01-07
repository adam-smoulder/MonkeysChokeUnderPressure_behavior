%%% This script makes the "main effect" figure from the behavior data. This
%%% is Figure 1C for now.
%%%
%%% Adam Smoulder, 7/23/20 (edit 7/27/20)

% Load subject data behavior files
loadBehaviorDataForMainFigures

%% Get trials by direction, reward, and trial status conditions
getNCondFromBehaviorForMainFigures
n_bySubject_byRew = cellfun(@(x) squeeze(sum(x,3))', cellfun(@(x) sum(x), n_byCond_bySubject,'uniformoutput',false),'uniformoutput',false);
n_bySubject_byRew = [n_bySubject_byRew{:}]'; % put in matrix form

%% Do bootstrapping and binomial proportion test based on these counts
nboots = 10000;
succRate_mean_bySubject_byReward = nan(nsubjects,nrewards);
succRate_bar95low_bySubject_byReward = nan(nsubjects,nrewards);
succRate_bar95up_bySubject_byReward = nan(nsubjects,nrewards);
succRate_sem_bySubject_byReward = nan(nsubjects,nrewards);
pvals_byRewardXReward_bySubject = nan(nrewards,nrewards,nsubjects);
for f = 1:nsubjects
    curCounts = n_byCond_bySubject{f};
    eventCounts = squeeze(sum(curCounts(:,:,1),1)); % successes are the first page (dim 3); we marginalize over directions
    totalCounts = squeeze(sum(sum(curCounts,3),1)); % all trials by reward size
    [~, bootMean, ~, ~, bar95low, bar95up, bootSem] = ...
        bootstrapBinaryEvent(eventCounts,totalCounts,nboots);
    succRate_mean_bySubject_byReward(f,:) = bootMean./totalCounts'*100;
    succRate_bar95low_bySubject_byReward(f,:) = bar95low./totalCounts'*100;
    succRate_bar95up_bySubject_byReward(f,:) = bar95up./totalCounts'*100;
    succRate_sem_bySubject_byReward(f,:) = bootSem./totalCounts'*100;
    
    % Use the binomial proportion test for significance; not related to
    % bootstrap, but might as well nab it here
    pvals_byRewardXReward_bySubject(:,:,f) = binomialProportionTest(eventCounts,totalCounts);
    
    disp(['Completed bootstrapping success rates for ' subjectNames{f}])
end; clear f




%% Make the plot for success rates with the SEM
rewNames = {'Small','Medium','Large','Jackpot'};

figure
for f = 1:nsubjects
    subplot(1,nsubjects,f)
    errorbar(rewards,succRate_mean_bySubject_byReward(f,:),...
       succRate_sem_bySubject_byReward(f,:),...
       'k.-','linewidth',0.5,'markersize',20)
   axis([0.5 0.5+nrewards 37.5 92.5])
   set(gca,'fontname','arial')
   set(gca,'fontsize',16)
   xticks(1:nrewards)
   xticklabels(rewNames)
   xtickangle(45)
   yticks([45 65 85])
   if f == 1
       ylabel('Success Rate (%)')
   end
   title(['Monkey ' subjectNames{f}(1) ', ' num2str(length(unique([behavior_bySubject{f}.day]))) ' sessions'])
end; clear f
set(gcf,'position',[418 548 1042 400])


pvals_byRewardXReward_bySubject


% Save it!
figname = 'Fig1B_mainEffect';
saveas(gcf,['MATLABFigs\' figname])
saveas(gcf,['SVGs\' figname '.svg'])


%% display p values and differences for S->L and L->J success rates
n_bySubject_byRew

for f = 1:nsubjects
    disp(['Subject ' subjectNames{f}(1)])
        disp(['S->L p = ' num2str(pvals_byRewardXReward_bySubject(1,3,f))])
        disp(['L->J p = ' num2str(pvals_byRewardXReward_bySubject(3,4,f))])
end; clear f

disp(diff(succRate_mean_bySubject_byReward(:,[1 3 4]),[],2))
