%%% This script adds the structure "reach_timings" to behavior. This
%%% contains 4 values: reactionTime (based on exit from center target),
%%% earlyTime (based on 1/3 of reach distance), middleTime (based on 1/3 -
%%% 2/3 of reach distance), and homingTime (based on 2/3 of reach distance
%%% to 1mm away from target). 
%%%
%%% Use the first line to specify the file to update.
%%%
%%% Adam Smoulder, 9/14/20

% fname = 'D:\AdamMatlab\~chokingBehaviorManuscript\Nelson\Nelson_behaviorFeatures_24-Jul-2020092033_r9';
% fname = 'D:\AdamMatlab\~chokingBehaviorManuscript\Ford\No Tubes\Ford_behaviorFeatures_24-Jul-2020092528_r8';
fname = 'D:\AdamMatlab\~chokingBehaviorManuscript\Earl_extendedRewardRange_forPaper\EarlExtRew_behaviorFeatures_07-Sep-2020134023_r8';
load(fname)
homingStopDistance = 1; % mm - how far from end target to stop homing time calculation
homingStopAddedTime = 150; % ms - how many ms beyond fail time to allow for homing time calc


%% Run the thing

startTargRad = behavior(1).task_parameters.startTargRadius;
reachTargRad = behavior(1).task_parameters.endTargRadius;
failTime = behavior(1).task_parameters.failTime;
targetDistance = behavior(1).task_parameters.targetDistance;
middleStartDistance = (targetDistance-startTargRad-reachTargRad)*2/3+reachTargRad % when we get within 85mm of end
homingStartDistance = (targetDistance-startTargRad-reachTargRad)*1/3+reachTargRad % use 1/3 left to target
homingEndDistance = homingStopDistance+reachTargRad;
middleEndDistance = homingStartDistance;
earlyEndDistance = middleStartDistance;
homingMaxEndTime = failTime+homingStopAddedTime;
for i = 1:length(behavior)
    curBeh = behavior(i);
    behavior(i).reach_timings.reactionTime = nan;
    behavior(i).reach_timings.earlyTime = nan;
    behavior(i).reach_timings.middleTime = nan;
    behavior(i).reach_timings.homingTime = nan;
    if ismember(curBeh.trialStatusLabels,[1 -22 -23 -24 -25 -31 -32 -34 -35]) % no reach times matter if not these listed classes
        
        % Get reaction time
        goCueTime = curBeh.t3_goCueTimes;
        kinTime = curBeh.kinematics_updated.time-goCueTime;
        distsFromEnd = curBeh.kinematics_updated.distanceFromEndTarget;
        distAlongTargetAxis = curBeh.kinematics_updated.rotatedPosition(:,1);
        behavior(i).reach_timings.reactionTime = curBeh.t4c_reactionTimes_exit-goCueTime;
        if behavior(i).reach_timings.reactionTime < 50 || behavior(i).reach_timings.reactionTime > failTime-100 % too quick/slow of reaction or bug, shouldn't be here. I do fail time -100 bc the fastest success reaches the target ~140ms after RT, so conservatively, we need at least 100ms
            disp(['Bad RT, subject goodBeh ind ' num2str(i) ' RT = ' num2str(behavior(i).reach_timings.reactionTime)])
            behavior(i).reach_timings.reactionTime = nan;
            continue
        end
        timeSinceGC = behavior(i).reach_timings.reactionTime;
        
        
        % Get early time
        upperEarlyInd = find((distsFromEnd <= earlyEndDistance) & (kinTime >= timeSinceGC),1);
        lowerEarlyInd = upperEarlyInd-1;
        if isempty(upperEarlyInd)
            disp(['No early completed, subject goodBeh ind ' num2str(i)])
            continue
        else
            behavior(i).reach_timings.earlyTime = twoPointInterpolation(earlyEndDistance,...
                distsFromEnd(lowerEarlyInd),distsFromEnd(upperEarlyInd),...
                kinTime(lowerEarlyInd),kinTime(upperEarlyInd))...
                -timeSinceGC;
        end
        timeSinceGC = timeSinceGC+behavior(i).reach_timings.earlyTime;
        
        % If early time is beyond trial limits, exclude
        if timeSinceGC > failTime
            disp(['Super late early completion, goodBeh ind ' num2str(i), ', ET+RT =  ' num2str(behavior(i).reach_timings.earlyTime+behavior(i).reach_timings.reactionTime)])
            behavior(i).reach_timings.earlyTime = nan;
            continue
        end
        
        
        % Get middle time
        upperMiddleInd = find((distsFromEnd <= middleEndDistance) & (kinTime >= timeSinceGC),1);
        lowerMiddleInd = upperMiddleInd-1;
        if isempty(upperMiddleInd)
            disp(['No middle completed, goodBeh ind ' num2str(i)])
            continue
        else
            behavior(i).reach_timings.middleTime = twoPointInterpolation(earlyEndDistance,...
                distsFromEnd(lowerMiddleInd),distsFromEnd(upperMiddleInd),...
                kinTime(lowerMiddleInd),kinTime(upperMiddleInd))...
                -timeSinceGC;
        end
        timeSinceGC = timeSinceGC+behavior(i).reach_timings.middleTime;
        
        % If middle time is beyond trial limits, exclude
        if behavior(i).reach_timings.middleTime+behavior(i).reach_timings.earlyTime+behavior(i).reach_timings.reactionTime > failTime
            disp(['Super late middle completion, goodBeh ind ' num2str(i), ', MT+ET+RT =  ' num2str(behavior(i).reach_timings.middleTime+behavior(i).reach_timings.earlyTime+behavior(i).reach_timings.reactionTime)])
            behavior(i).reach_timings.middleTime = nan;
            continue
        end
        
        
        
        % Get homing time
        upperHomingInd = find((distsFromEnd <= homingEndDistance) & (kinTime >= timeSinceGC),1);
        %             if isempty(upperHomingInd) % may be wide-overshoot, so never got close enough; call it when they pass by the midway of the target
        %                 upperHomingInd = find(distAlongTargetAxis >= targetDistance,1);
        %             end
        lowerHomingInd = upperHomingInd-1;
        if isempty(upperHomingInd)
            disp(['No homing achieved, goodBeh ind ' num2str(i)])
            continue
        else
            behavior(i).reach_timings.homingTime = twoPointInterpolation(homingEndDistance,...
                distsFromEnd(lowerHomingInd),distsFromEnd(upperHomingInd),...
                kinTime(lowerHomingInd),kinTime(upperHomingInd))...
                -timeSinceGC;
        end
        timeSinceGC = timeSinceGC+behavior(i).reach_timings.homingTime;
        
        % If homing time is beyond trial limits + reaction
        % (homingStopAddedTime), exclude
        if timeSinceGC > homingMaxEndTime
            disp(['Homing too late, goodBeh ind ' num2str(i), ', HT+MT+ET+RT =  ' num2str(behavior(i).reach_timings.homingTime+behavior(i).reach_timings.middleTime+behavior(i).reach_timings.earlyTime+behavior(i).reach_timings.reactionTime)])
            behavior(i).reach_timings.homingTime = nan;
            continue
        end
        
    end
end; clear i


%% Overwrite
save(fname,'behavior')