%%% This script shows overshoots, undershoots, and reach epoch success 
%%% endpoints for Monkey N. For visualization/intuition, we remove some
%%% over/undershoots (see SI methods). See FigX_MonkeyNEndpoints_all for
%%% the version with all endpoints shown.
%%%
%%% Written by Patrick Marino, 10/18/20, edited by Adam Smoulder, 10/20/20

%% Load Data; Select Monkey N; Create Status Codes Vec
% Load subject data behavior files
loadBehaviorDataForMainFigures

% Get trials by direction, reward, and trial status conditions
getNCondFromBehaviorForMainFigures

% Only keep trials that make it to the reach; for reach epoch, we call
% anything that makes it to target hold (w/o overshooting) success (1), and
% have overshoots (-22) and undershoots (-23)
f = 1; % only using Nelson for this script, no need to loop
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


% Establish reward and status variables
rewards = 1:4;
rewardNames = {'Small','Medium','Large','Jackpot'};
numSubjects = length(subjectNames);
numRewards = length(rewardNames);

nelsonData = behavior_bySubject{1,1};
statusCodes = [1,-23,-22]; %Suc, Under, Over
numStatus = length(statusCodes);

%% Get Target Plotting Info
task_parameters = nelsonData(1).task_parameters;
r_end = task_parameters.endTargRadius;
endtarg = [task_parameters.targetDistance 0];
theta = linspace(0,2*pi,100);

%% Plot Monkey N Reach Endpoints For Direction 1 (right)
statusColors = [0 1 0; 1 1 1 ; 0 0 0 ] + [0 ; -0.3 ; 0.2];
dColor = [0 0 1]; %Diamond Color 
ms = 12; %Marker size
ds = 6; %Diamond Size
fs = 14; %Font size
f=figure('Position',[100 100 400 400]); hold on
f.Renderer = 'Painters';

directionInd = 1;
directionData = nelsonData([nelsonData.direction]==directionInd);
numJackpotTrials = sum([directionData.reward]==4);
jackpotEndpoints = NaN(numJackpotTrials,2);
jackpotTrialInd = 1;

% Plot reach target
plot(r_end*cos(theta)+endtarg(1),r_end*sin(theta)+endtarg(2),'k-')

%Plot endpoints according to status-specific rules; Save jackpot endpoints
%   to enclose with diamonds later

% We only want to plot N of each status; specify these as counters
counts = fliplr(histcounts([directionData.trialStatusLabels]));
counts(counts==0) = [];
countsper100 = round(counts./sum(counts)*100);
nhundred = 5; % how many hundreds of points to plot
nsucctoplot = countsper100(1)*nhundred;
nushoottoplot = countsper100(3)*nhundred;
noshoottoplot = countsper100(2)*nhundred;

succplotcount = 0;
ushootplotcount = 0;
oshootplotcount = 0;

% Randomize trial order
rng(3195) % seed rng for repeatability

for status = statusCodes
    directionStatusData = directionData([directionData.trialStatusLabels]==status); 
    numTrials = size(directionStatusData,2);
    shuffledTrials = randperm(numTrials);
    for shuffledTrialInd = 1:numTrials
        trial = shuffledTrials(shuffledTrialInd);
        
        %Get trajectory
        curKin = directionStatusData(trial).kinematics_updated;
        time = curKin.time;
        traj = curKin.rotatedPosition;
        reward = directionStatusData(trial).reward;
        
        %Get index of sample closest to reach end time
        reachEndTime = directionStatusData(trial).t6_reachEndTimes;
        relTime = time-reachEndTime;
        [~,reachEndInd] = min(abs(relTime));
        
        %Get index of sample closest to post reach state + 100ms
        postReachStateTime100 = directionStatusData(trial).t7_postReachStateTimes + 100;
        relTime = time-postReachStateTime100;
        [~,postReachStateInd100] = min(abs(relTime));
        
        %Get index of sample closest to trial end time
        trialEndTime = directionStatusData(trial).t8_trialEndTimes;
        relTime = time-trialEndTime;
        [~,trialEndInd] = min(abs(relTime));
        
        %Get index of sample closest to trial end time + 100ms
        trialEndTime100 = directionStatusData(trial).t8_trialEndTimes + 100;
        relTime = time-trialEndTime100;
        [~,trialEndInd100] = min(abs(relTime));
        
        %Plot
        switch status
            %For success, plot endpoints that fell inside the target 
            case 1
                xCoord = traj(reachEndInd,1); yCoord = traj(reachEndInd,2);
                %Uncomment next line to instead plot endpoints 100ms after target entry
                %   xCoord = traj(postReachStateInd100,1); yCoord = traj(postReachStateInd100,2);
                r = sqrt((xCoord-endtarg(1))^2 + (yCoord-endtarg(2))^2);
                if r < r_end && succplotcount < nsucctoplot
                    plot(xCoord,yCoord,'.','MarkerSize',ms,'Color',statusColors(1,:));
                    %If trial was a jackpot, save the endpoint location 
                    if reward == 4
                        jackpotEndpoints(jackpotTrialInd,:) = [xCoord,yCoord];
                        jackpotTrialInd = jackpotTrialInd + 1;
                    end
                    succplotcount = succplotcount+1;
                end
                
            %For undershoot, plot endpoints that fell outside of target and left of
            %target center
            case -23
                xCoord = traj(reachEndInd,1); yCoord = traj(reachEndInd,2);
                %Uncomment next line to instead plot endpoints at fail time
                %   xCoord = traj(trialEndInd,1); yCoord = traj(trialEndInd,2);
                r = sqrt((xCoord-endtarg(1))^2 + (yCoord-endtarg(2))^2);
                if r > r_end && xCoord < endtarg(1) && ushootplotcount < nushoottoplot
                    plot(xCoord,yCoord,'.','MarkerSize',ms,'Color',statusColors(2,:));
                    %If trial was a jackpot, save the endpoint location 
                    if reward == 4
                        jackpotEndpoints(jackpotTrialInd,:) = [xCoord,yCoord];
                        jackpotTrialInd = jackpotTrialInd + 1;
                    end
                    ushootplotcount = ushootplotcount+1;
                end
                
            %For overshoot, plot endpoints that fell outside of target and right of 
            %target center; technically, this is NOT how all overshoots
            %necessarily have to go (i.e., you could go through the target
            %at some angle such that you don't end up necessarily on the
            %opposite side of where the center target is), but this
            %figure's goal is to primarily serve as schematic.
            case -22
                xCoord = traj(reachEndInd,1); yCoord = traj(reachEndInd,2);
                %Uncomment next line to instead plot endpoints 100ms after fail time
                %   xCoord = traj(trialEndInd100,1); yCoord = traj(trialEndInd100,2);
                r = sqrt((xCoord-endtarg(1))^2 + (yCoord-endtarg(2))^2);
                if r > r_end && xCoord > endtarg(1) && oshootplotcount < noshoottoplot
                    plot(xCoord,yCoord,'.','MarkerSize',ms,'Color',statusColors(3,:));
                    %If trial was a jackpot, save the endpoint location 
                    if reward == 4
                        jackpotEndpoints(jackpotTrialInd,:) = [xCoord,yCoord];
                        jackpotTrialInd = jackpotTrialInd + 1;
                    end
                    oshootplotcount = oshootplotcount+1;
                end
        end
    end
end

%Not all jackpot trials met the rules for plotting.  So get rid of 
%   placeholders for those which didn't
jackpotEndpoints = jackpotEndpoints(~isnan(jackpotEndpoints(:,1)),:);
numJackpotTrials = size(jackpotEndpoints,1);

%Plot diamonds around jackpot endpoints 
for trial = 1:numJackpotTrials
   xCoord = jackpotEndpoints(trial,1);
   yCoord = jackpotEndpoints(trial,2);
   plot(xCoord,yCoord,'d','MarkerSize',ds,'Color',dColor); 
end



%Set up plot axes, ticks, labels, and title
axis equal
axis([72.5 95 -10 10])
set(gca,'fontname','arial')
set(gca,'fontsize',fs)
set(gca,'TickDir','out')
xlabel('On-Target Axis (mm)')
ylabel('Off-Target Axis (mm)')
xticks([75 95])
xticklabels({75 95})
yticks([-10 10])
yticklabels({-10 10})
title(['Monkey N ' num2str(100*nhundred) ' Reach Endpoints'])

% Save it!
figname = 'Fig4A_monkeyNEndpoints';
saveas(gcf,['MATLABFigs\' figname])
saveas(gcf,['SVGs\' figname '.svg'])

size(jackpotEndpoints)