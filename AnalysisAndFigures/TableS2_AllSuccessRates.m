%%% This script compiles all the data from the individual days of the
%%% normal task and combined across sessions data from controls/other
%%% tasks into a format easily copyable to make Table S2 in Excel
%%%
%%% Adam Smoulder, 8/10/20 (edit 10/20/20)

% Load subject data behavior files
subjectNames = {...
    'Monkey N',...
    'Monkey F',...
    'Monkey E',...
    'Monkey E standard reward range',...
    'Monkey F precision Task',...
    'Monkey F precision task, Mini-Jackpot sessions',...
    'Monkey E Rare-Large sessions',...
    'Monkey E Common-Jackpot sessions',...
    };
behaviorFnames = {...
    'DataFolder\Nelson_behaviorFeatures_20210107',...
    'DataFolder\Ford_behaviorFeatures_20210107',...
    'DataFolder\Earl_behaviorFeatures_20210107',...
    'DataFolder\EarlStandardRewRange_behaviorFeatures_20210107',...
    'DataFolder\FordTubes_minimumTubeFeatures_27-Jul-2020133051',...
    'DataFolder\FordTubesMJ_minimumTubeFeatures_27-Jul-2020133443',...
    'DataFolder\EarlRL_behaviorFeatures_20210107',...
    'DataFolder\EarlCJ_behaviorFeatures_20210107',...
    };
nsubjects = length(subjectNames);
behavior_bySubject = cell(nsubjects,1);
for f = 1:nsubjects
    load(behaviorFnames{f},'behavior');
    behavior_bySubject{f} = behavior;
    disp(['Loaded subject ' subjectNames{f}])
end; clear f
clear behavior


%% For each subject, we count total number of valid (non-quitout) trials
%  and successes; for the first 3 subjects (main data), we also split this
%  by day. For the others, we call it effectively 1 day.

% For bad statuses
%  - quitouts and misstarts (no real trial) are excluded ([-20 -13 -14 0] = 0)
%  - No attempts, wild, and early returns are also removed due to infrequency ([-21 -26 -33] = 0)
%  - the few -111s we have are misstarts (-111 = 0)
%  - we'll remove any excluded trials (0s removed)
% So for any nontubes task, [0 -111 -13 -14 -20 -21 -26 -33] are bad

n_bySubject_byRew_byDay = cell(nsubjects,1);
nsucc_bySubject_byRew_byDay = cell(nsubjects,1);
for f = 1:nsubjects
    % Precision task data has 0 as failure, while others have 0 as quitout,
    % so differentiate:
    if ismember(f,[5,6])
        badStatuses = [];
    else
        badStatuses = [0 -111 -13 -14 -20 -21 -26 -33];
    end
    
    % Get the current behavior and relevant labels
    curBeh = behavior_bySubject{f};
    trialStatusLabels = [curBeh.trialStatusLabels];
    rewardLabels = [curBeh.reward];
    rewards = unique(rewardLabels);
    nrewards = length(rewards);
    if f <= 3 % Main data - show all days
        dayLabels = [curBeh.day];
        days = unique(dayLabels);
        ndays = length(days);
    else
        dayLabels = ones(size(trialStatusLabels)); % just call it all 1 day
        days = 1;
        ndays = 1;
    end
    
    % and do the counting
    n_bySubject_byRew_byDay{f} = nan(nrewards,ndays);
    nsucc_bySubject_byRew_byDay{f} = nan(nrewards,ndays);
    for r = 1:nrewards
        for a = 1:ndays
            n_bySubject_byRew_byDay{f}(r,a) = sum(...
                      rewardLabels==rewards(r) & ...
                      dayLabels==days(a) & ...
                      ~ismember(trialStatusLabels,badStatuses));
            nsucc_bySubject_byRew_byDay{f}(r,a) = sum(...
                      rewardLabels==rewards(r) & ...
                      dayLabels==days(a) & ...
                      trialStatusLabels==1);
        end; clear a
    end; clear r
    
end; clear f

%% Now we make a big fat table...to do this, we need to fill out row 
%  info for each of the columns. The columns are going to be:
%  - Subject Name (name)
%  - Day number (dayNum)
%  - Small Reward Successful/Total Trials (SRR)
%  - Small Reward Stats (SRS)
%  - Medium Reward Successful/Total Trials (MRR)
%  - Medium Reward Stats (MRS)
%  - Large Reward Successful/Total Trials (LRR)
%  - Large Reward Stats (LRS)
%  - Jackpot Reward Successful/Total Trials (JRR)
%  - Jackpot Reward Stats (JRS)
%  - Rare-Large/Common-Jackpot Reward Successful/Total Trials (ERR)
%  - Rare-Large/Common-Jackpot Reward Stats (ERS)
%
% Each day gets 4 rows. In this, the subject name and day are repeated, the
% top 2 rows have the % success while bottom 2 have #succ/#total, and each
% row is the sig stars vs the other reward sizes. The exception is for the
% RL/CJ data, which has 5 rows (bottom 3 are #)
%
% ALL of the data in this is text.

bigTableData = cell(12,1); % one for each column
maxNrewards = 5; % highest number of rewards

for f = 1:nsubjects
    curN_byRew_byDay = n_bySubject_byRew_byDay{f};
    curNsucc_byRew_byDay = nsucc_bySubject_byRew_byDay{f};
    [nrewards,ndays] = size(curN_byRew_byDay);
    curName = subjectNames{f}
    
    for a = 1:ndays
        curN_byRew = curN_byRew_byDay(:,a);
        curNSucc_byRew = curNsucc_byRew_byDay(:,a);
        curRates = curNSucc_byRew./curN_byRew;
        dayPvalMatrix = binomialProportionTest(curNSucc_byRew,curN_byRew);
        
        for r1 = 1:nrewards     % This will be for rows
            bigTableData{1}(end+1) = {curName};             % Column 1: Subject Name
            bigTableData{2}(end+1) = {['Day ' num2str(a)]}; % Column 2: Day
            
            for r2 = 1:maxNrewards  % this will be for columns
                if r2 > nrewards % no rewards left; fill with blanks
                    bigTableData{1+2*r2}(end+1) = {''};
                    bigTableData{2+2*r2}(end+1) = {''};
                    continue
                end
                
                % Column 3 5 7 9 11: Successful/Total Trials 
                if r1 < 3 % top 2 rows of this reward show percentage success
                    bigTableData{1+2*r2}(end+1) = ...
                        {[num2str(round(curNSucc_byRew(r2)/curN_byRew(r2)*100,1)) '%']};
                else % below that shows the #succ / # fail
                    bigTableData{1+2*r2}(end+1) = ...
                        {[num2str(curNSucc_byRew(r2)) '/' num2str(curN_byRew(r2))]};
                end
                
                % Column 4 6 8 10 12: Stats
                if r1 == r2 % use dashes for stats
                    bigTableData{2+2*r1}(end+1) = {'-'};
                else % use stars based on p value
                    curPVal = dayPvalMatrix(r1,r2);
                    if curPVal < 0.001
                        sigStr = '***';
                    elseif curPVal < 0.01
                        sigStr = '**';
                    elseif curPVal < 0.05
                        sigStr = '*';
                    else
                        sigStr = 'ns';
                    end
                    bigTableData{2+2*r2}(end+1) = {sigStr};
                    if f ==4
                        curPVal
                        disp('')
                    end
                end
            end; clear r2
        end; clear r1
    end; clear a
end; clear f

% And finally make it into a real table and export!
bigTable = table(...
    bigTableData{1}',...
    bigTableData{2}',...
    bigTableData{3}',...
    bigTableData{4}',...
    bigTableData{5}',...
    bigTableData{6}',...
    bigTableData{7}',...
    bigTableData{8}',...
    bigTableData{9}',...
    bigTableData{10}',...
    bigTableData{11}',...
    bigTableData{12}'...
    );

tableName = 'TableS2_AllSuccessRates';
writetable(bigTable,[tableName '.csv'])


%  NOTE the output table for RL and CJ stuff has them before Jackpot, so
%  the columns won't line up exactly; you have to flip-flop those
