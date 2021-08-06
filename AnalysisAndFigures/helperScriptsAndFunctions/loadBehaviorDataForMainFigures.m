%% This script loads the relevant behavioral files for each subject for the
%  main figures of the paper.
%
% Adam Smoulder, 7/23/20 (edit 10/15/20)

% Load subject data behavior files
subjectNames = {'Nelson','Ford','Earl'};
behaviorFnames = {...
    'Nelson_behaviorFeatures_20210107',...
    'Ford_behaviorFeatures_20210107',...
    'Earl_behaviorFeatures_20210107',...
    };
nsubjects = length(subjectNames);
behavior_bySubject = cell(nsubjects,1);
for f = 1:nsubjects
    load(behaviorFnames{f},'behavior')
    behavior_bySubject{f} = behavior;
    disp(['Loaded subject ' subjectNames{f}])
end; clear f
clear behavior