%%% This script looks at success rates as a function of day x target. We'll
%%% just use Monkey E for this example, since he has the most data.
%%%
%%% Adam Smoulder, 11/10/20

% Load subject data behavior files
loadBehaviorDataForMainFigures

%% Get trials by direction, reward, and trial status conditions
getNCondFromBehaviorForMainFigures

%% We need to re-count and split by reward x direction x day. We
%  only care about success versus failure, to simplify things
nsplits = 4;

for f = 1:nsubjects
    behavior = behavior_bySubject{f};
    trialStatusLabels = trialStatusLabels_bySubject{f};
    rewardLabels = rewardLabels_bySubject{f};
    dirLabels = dirLabels_bySubject{f};
    dirs = unique(dirLabels);
    ndirs = length(dirs);
    dayLabels = [behavior.day]';
    days = unique(dayLabels);
    ndays = length(days);
    
    n_byRew_byDir_byDay = nan(nrewards,ndirs,ndays);
    nsucc_byRew_byDir_byDay = nan(nrewards,ndirs,ndays);
    
    for d = 1:ndirs
        for r = 1:nrewards
            for a = 1:ndays
                curInds = dirLabels==dirs(d) &...
                    rewardLabels==rewards(r) &...
                    dayLabels==days(a);
                n_byRew_byDir_byDay(r,d,a) = sum(curInds);
                nsucc_byRew_byDir_byDay(r,d,a) = sum(curInds & trialStatusLabels==1);
            end; clear a
        end; clear r
    end; clear d
    
    % And combine days with dirs
    n_byRew_byDayDir = n_byRew_byDir_byDay(:,:);
    nsucc_byRew_byDayDir = nsucc_byRew_byDir_byDay(:,:);
    
    
    % First, order by median success rates, then split into quartiles
    [~,order] = sort(mean(nsucc_byRew_byDayDir([2:3],:)./n_byRew_byDayDir([2:3],:),1),'descend');
    
    n_byRew_byPrctile = nan(nrewards,nsplits);
    nsucc_byRew_byPrctile = nan(nrewards,nsplits);
    count = 0;
    for n = 1:nsplits
        curInds = (count+1):(count+round(ndirs*ndays/nsplits));
        n_byRew_byPrctile(:,n) = sum(n_byRew_byDayDir(:,order(curInds)),2);
        nsucc_byRew_byPrctile(:,n) = sum(nsucc_byRew_byDayDir(:,order(curInds)),2);
        count = count+round(ndirs*ndays/nsplits);
    end; clear n
    
    nsucc_byRew_byPrctile'./n_byRew_byPrctile'
    
    n_bySubject_byRew_byPrctile(f,:,:) = n_byRew_byPrctile;
    nsucc_bySubject_byRew_byPrctile(f,:,:) = nsucc_byRew_byPrctile;

end; clear f

%% We want to compare the top and bottom quartiles mainly - so let's just
%  combine the middle 50% and have it there
if nsplits > 3
    n_bySubject_byRew_byPrctile(:,:,2) = sum(n_bySubject_byRew_byPrctile(:,:,2:end-1),3); % sum all of the middle ones
    n_bySubject_byRew_byPrctile(:,:,3) = n_bySubject_byRew_byPrctile(:,:,end);            % keep the last prctile
    n_bySubject_byRew_byPrctile(:,:,4:end) = []; %  remove all else.
    nsucc_bySubject_byRew_byPrctile(:,:,2) = sum(nsucc_bySubject_byRew_byPrctile(:,:,2:end-1),3); % sum all of the middle ones
    nsucc_bySubject_byRew_byPrctile(:,:,3) = nsucc_bySubject_byRew_byPrctile(:,:,end);            % keep the last prctile
    nsucc_bySubject_byRew_byPrctile(:,:,4:end) = []; %  remove all else.
    nsplits = 3;
end

%% Bootstrap
nboots = 10000;
succRate_mean_bySubject_byReward_byPrctile = nan(nsubjects,nrewards,nsplits);
succRate_bar95low_bySubject_byReward_byPrctile = nan(nsubjects,nrewards,nsplits);
succRate_bar95up_bySubject_byReward_byPrctile = nan(nsubjects,nrewards,nsplits);
succRate_sem_bySubject_byReward_byPrctile = nan(nsubjects,nrewards,nsplits);
pvals_byRewardXReward_bySubject_byPrctile = nan(nrewards,nrewards,nsubjects,nsplits);
for f = 1:nsubjects
    for n = 1:nsplits
        eventCounts = nsucc_bySubject_byRew_byPrctile(f,:,n); % successes are the first page (dim 3); we marginalize over directions
        totalCounts = n_bySubject_byRew_byPrctile(f,:,n); % all trials by reward size
        [~, bootMean, ~, ~, bar95low, bar95up, bootSem] = ...
            bootstrapBinaryEvent(eventCounts,totalCounts,nboots);
        succRate_mean_bySubject_byReward_byPrctile(f,:,n) = bootMean./totalCounts'*100;
        succRate_bar95low_bySubject_byReward_byPrctile(f,:,n) = bar95low./totalCounts'*100;
        succRate_bar95up_bySubject_byReward_byPrctile(f,:,n) = bar95up./totalCounts'*100;
        succRate_sem_bySubject_byReward_byPrctile(f,:,n) = bootSem./totalCounts'*100;
        
        % Use the binomial proportion test for significance; not related to
        % bootstrap, but might as well nab it here
        pvals_byRewardXReward_bySubject_byPrctile(:,:,f,n) = binomialProportionTest(eventCounts,totalCounts);
    end; clear n
    disp(['Completed bootstrapping success rates for ' subjectNames{f}])
end; clear f


%% 
rewNames = {'Small','Medium','Large','Jackpot'};
figure
for f = 1:nsubjects
    subplot(1,3,f); 
    hold on
    topMeans = succRate_mean_bySubject_byReward_byPrctile(f,:,1);
    topSems = succRate_sem_bySubject_byReward_byPrctile(f,:,1);
    middleMeans = succRate_mean_bySubject_byReward_byPrctile(f,:,ceil(nsplits/2));
    middleSems = succRate_sem_bySubject_byReward_byPrctile(f,:,ceil(nsplits/2));
    bottomMeans = succRate_mean_bySubject_byReward_byPrctile(f,:,end);
    bottomSems = succRate_sem_bySubject_byReward_byPrctile(f,:,end);
    errorbar(rewards-0.15,topMeans,bottomSems,'r.-','linewidth',1,'markersize',20)
    errorbar(rewards+0.15,bottomMeans,bottomSems,'b.-','linewidth',1,'markersize',20)
    errorbar(rewards,middleMeans,middleSems,'-','linewidth',1,'markersize',20,'color',0.6*ones(1,3),'capsize',1)
    if f == 1
        legend('Top 25%','Bottom 25%','Middle 50%','location','SW')
    end
    axis([0.5 4.25 27.5 95])
    yticks([30 60 90])
    if f == 1
       ylabel('Success Rate (%)')
    end
   set(gca,'fontname','arial')
   set(gca,'fontsize',16)
   xticks(1:nrewards)
   xticklabels(rewNames)
   xtickangle(45)
   title(['Monkey ' subjectNames{f}(1) ', ' num2str(length(unique([behavior_bySubject{f}.day]))) ' sessions'])
end; clear f
set(gcf,'position',[418 548 1042 400])


% save it!
figname = 'FigS7_difficultyByDayTarget';
saveas(gcf,['MATLABFigs\' figname])
saveas(gcf,['SVGs\' figname '.svg'])


