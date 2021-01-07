function [N_byCond, dirLabels, rewardLabels, trialStatusLabels, dirs,...
    rewards, statuses] = getNByCondFromBehavior(behavior)
% This function extracts the number of trials for each condition from a
% given behavior structure. When I say "condition" I'm referring to
% direction x reward x trial status (success or failure code).
%
% Input:
% - behavior: [ntrials x 1] structure with the fields "direction",
% "reward", and "trialStatusLabel" at minimum.
%
% Outputs:
% - N_byCond: [ndirections x nrewards x nstatuses] matrix with the counts
% for each condition.
% - dirLabels: [ntrials x 1] direction labels for each trial
% - rewardLabels: [ntrials x 1] reward labels for each trial
% - trialStatusLabels: [ntrials x 1] success or failure code for each trial
% - dirs: [ndirections x 1], unique(dirLabels)
% - rewards: [nrewards x 1], unique(rewardLabels)
% - statuses: [nstatuses x 1], flip(unique(trialStatusLabels)); we do the
% flip so that "1" (success) is first, followed by failures
%
% Adam Smoulder, 7/24/20
%
% CODES FOR EACH FAILURE MODE ARE AS FOLLOWS
%
% 1    % Success!
% 
% Delay epoch failures
% -11  % False start: attempted reach towards target early
% -12  % Cheat/drift: no full reach attempt made in any dir
% -13  % Mis-start: accidentally enters/exits start target
% -14  % Quit out: reach made away from target
% -111 % Shouldn't happen if it does, need to readjust criteria
% 
% Reach epoch failures
% -20 % Quit-out: Intentional reach in wrong direction
% -21 % No Attempt: no attempt was made (or a very lazy one, indicated by a super late peak speed even though the target was technically exited)
% -22 % Overshoot: reach missed the target
% -23 % Undershoot: main reach finished outside of the target
% -24 % Slow: Reach is too slow to have possibly made it
% -25 % Slow: The remainder are on the ascent/descent of the peak speed curve and are mid-reach.
% -26 % Wild: Reach is very very not straight!
% 
% Target hold epoch failures
% -31 % Scuff: Cursor just barely scrapes the target
% -32 % Overshoot: Reach goes straight through target
% -33 % Early Return: Reaching back to center before time is up
% -34 % Drift: Doesn't blow through target, but doesn't hold properly, drifting out of the target
% -35 % Jitter: Holding at edge of target -> jitters out

% Get all of the labels out that we need
dirLabels = [behavior.direction]';
dirs = unique(dirLabels);
ndirs = length(dirs);
rewardLabels = [behavior.reward]';
rewards = unique(rewardLabels);
nrewards = length(rewards);
trialStatusLabels = [behavior.trialStatusLabels]';
statuses = flip(unique(trialStatusLabels));
nstatuses = length(statuses);

% Run through each condition and sum the # of trials where it's true
N_byCond = nan(ndirs,nrewards,nstatuses);
for d = 1:ndirs
    for r = 1:nrewards
        for s = 1:nstatuses
            N_byCond(d,r,s) = sum(dirLabels==dirs(d) &...
                rewardLabels==rewards(r) &...
                trialStatusLabels==statuses(s));
        end; clear s
    end; clear r
end; clear d

end

