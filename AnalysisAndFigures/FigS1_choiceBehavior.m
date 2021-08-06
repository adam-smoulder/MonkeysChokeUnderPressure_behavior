%%% This script makes the figure to show choice behavior. While our
%%% original experiments only had left/right choices (which made for really
%%% easy/pretty plots), Earl's sessions have all 8, so we'll split into 2
%%% categories based on his direction of bias (see text, Fig S1A)
%%%
%%% Adam Smoulder, 8/11/20 (edit 8/21/20)

% Load choice behavior for each subject
subjectNames = {...
    'Monkey N',...
    'Monkey F',...
    'Monkey E',...
    'Monkey E Rare-Large',...
    'Monkey E Common-Jackpot',...
    };
behaviorFnames = {...
    'SAVESTATE_ChoiceBehaviorExt_Nelson',...
    'SAVESTATE_ChoiceBehaviorExt_Ford',...
    'SAVESTATE_ChoiceBehaviorExt_Earl',...
    'SAVESTATE_ChoiceBehaviorExt_EarlRL',...
    'SAVESTATE_ChoiceBehaviorExt_EarlCJ',...
    };
nsubjects = length(subjectNames);
choiceBehavior_bySubject = cell(nsubjects,1);
for f = 1:nsubjects
    load(behaviorFnames{f},'choiceBehavior');
    choiceBehavior_bySubject{f} = choiceBehavior;
    disp(['Loaded subject ' subjectNames{f}])
end; clear f
clear choiceBehavior

%% We need to populate Red (right) vs blue (left) choices. We'll count the
%  number of right choices. Red is directions 1-4, Blue is 5-8.

% We assume that rewards are 1:nrewards in terms of labels. This will break
% if not!
ntotal_bySubject_byRewR_byRewB = cell(nsubjects,1); % R = red, B = blue
nR_bySubject_byRewR_byRewB = cell(nsubjects,1);

% There's a lot of ways we could do this, but I'll just do the
% inefficient-ass way of going through every trial and putting it in the
% right matrix
for f = 1:nsubjects
    curBeh = choiceBehavior_bySubject{f};
    nrewards = length(curBeh.comboN);
    ntrials = length(curBeh.decisionChoices);
    ntotal_bySubject_byRewR_byRewB{f} = zeros(nrewards,nrewards);
    nR_bySubject_byRewR_byRewB{f} = zeros(nrewards,nrewards);
    
    for n = 1:ntrials
        rewOps = curBeh.rewardOptions(n,:);
        dirOps = curBeh.directionOptions(n,:);
        choice = curBeh.decisionChoices(n); % 0 for first choice, 1 for second choice (confusing, I know)
        
        % We're going to reorganize so red target is always first
        if dirOps(2) < dirOps(1) % out of order, flip all options
            dirOps = fliplr(dirOps);
            rewOps = fliplr(rewOps);
            choice = 1-choice;
        end
        
        % Add this trial to the appropriate counter
        ntotal_bySubject_byRewR_byRewB{f}(rewOps(1),rewOps(2)) = ...
            1+ntotal_bySubject_byRewR_byRewB{f}(rewOps(1),rewOps(2));
        if choice==0 % if red choice
            nR_bySubject_byRewR_byRewB{f}(rewOps(1),rewOps(2)) = ...
                1+nR_bySubject_byRewR_byRewB{f}(rewOps(1),rewOps(2));
        end
    end; clear n
end; clear f


%% Make figures
figure
set(gcf,'position',[-3803 63 3693 421])
% red = [236 34 36]/255; % 100% red choice
% blue = [67 120 187]/255; % 100% blue choice
red = [255 0 0]/255;
blue = [0 0 255]/255;
midpoint = [255 255 255]/255; % 50/50
r = [linspace(red(1),midpoint(1),32) linspace(midpoint(1),blue(1),32)]';
g = [linspace(red(2),midpoint(2),32) linspace(midpoint(2),blue(2),32)]';
b = [linspace(red(3),midpoint(3),32) linspace(midpoint(3),blue(3),32)]';
rbColormap = flip([r g b]);

% [236 34 36] ec2224
% to
% [67 120 187] 4377bb

for f = 1:nsubjects
    curBeh = choiceBehavior_bySubject{f};
    nrewards = length(curBeh.comboN);
    if f <= 3 % just 1-4
        rewNames = {'S','M','L','J'};
    elseif f == 4
        rewNames = {'S','M','L','RL','J'};
    elseif f == 5
        rewNames = {'S','M','L','CJ','J'};
    end
    
    subplot(1,nsubjects,f)
    vals = (nR_bySubject_byRewR_byRewB{f}./ntotal_bySubject_byRewR_byRewB{f})';
    vals(isnan(vals)) = 0.5; % just call it 50/50 if not shown for color's sake; text will show 0/0
    imagesc(1:nrewards,-1:-1:-nrewards,vals)
    axis xy
    hold on;
    c = colormap(rbColormap);
    xticks(1:nrewards)
    xticklabels(rewNames)
    yticks(-nrewards:-1)
    yticklabels(flip(rewNames))
    colorbar
    set(gca,'clim',[0 1])
    set(gca,'fontsize',12)
    set(gca,'fontname','arial')
    title(subjectNames{f})
    
    % Calculate how many of the incorrect choices were in biased directions
    biasDirs = find(curBeh.sameRewDirFracs > 0.5);
    
    nCorrect = nansum(curBeh.comboNCorrect(:));
    nWrongInBias = sum(curBeh.wrongDecSelectedDirs(biasDirs));
    nWrongOutBias = nansum(curBeh.wrongDecSelectedDirs(:))-nWrongInBias;
    
    disp(subjectNames{f})
    nR_bySubject_byRewR_byRewB{f}'
    ntotal_bySubject_byRewR_byRewB{f}'
end; clear f



% Save it!
figname = 'FigS1_choiceBehavior';
saveas(gcf,['MATLABFigs\' figname])
saveas(gcf,['SVGs\' figname '.svg'])



%% For the main text, get the total % correct overall
totalN = 0;
totalNCorrect = 0;
for f = 1:3
    curComboN = choiceBehavior_bySubject{f}.comboN;
    curComboNCorrect = choiceBehavior_bySubject{f}.comboNCorrect;
    
    % we have to skip the diagonal, since it tells us nothing about
    % "correct" choices
    totalN = totalN+nansum(curComboN-diag(diag(curComboN)));
    totalNCorrect = totalNCorrect+nansum(curComboNCorrect-diag(diag(curComboNCorrect)));
end; clear f

round(100*totalNCorrect/totalN,1)