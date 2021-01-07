%%% This is to make part of the schematic for the BREP/homing time figure.
%%% Mostly it summons the data for the example trial and plots some of the
%%% key points, and the rest of the schematic making is done in
%%% illustrator.
%%%
%%% Adam Smoulder, 8/7/20 (edit: 1/7/21)

%% Load stuff
% for schematic
load('SAVESTATE_forBREPExample2')

%% Example schematic BREP figure
figure
set(gcf,'position',[7 29 834 938])
startTime = 1304;
endTime = 1950;

GOLD = [1,0.647,0]*0.9;

% speed trace
subplot(4,2,[1 3])
firstInd = find(curTime<=startTime,1,'last');
lastInd = find(curTime>=endTime,1,'first');
indsToCompare = compareInds;
hold on
plot(curTime(1:lastIndToPlot)-startTime,speed(1:lastIndToPlot),'k-','linewidth',2)
plot(curTime(indsToCompare(7:end-7))-startTime,sqrt(sum(ballisticPredVel(7:end-7,:).^2,2)),'k--','linewidth',1.5,'color',GOLD)
plot(curTime(lastInd-7)-startTime,sqrt(sum(ballisticPredVel(end-7,:).^2,2)),'k+','linewidth',2,'markersize',10,'color',[0.65 0 0])
plot(371*ones(1,100), linspace(0,0.4274,100),'k-','color',GOLD,'linewidth',0.5)
axis([0 endTime-startTime 0 inf])
xlabel('Time since go cue (ms)')
xticks([0 200 400 600])
yticks([0 0.1 0.2 0.3 0.4])
ylabel('Speed (m/s)')
set(gca,'fontsize',10)

% distance from end target trace
subplot(4,2,[2 4])
hold on
plot(curTime(1:lastIndToPlot)-startTime,distFromEndTarg(1:lastIndToPlot),'k-','linewidth',2)
plot(549,sqrt((85-curRotatedBallisticPredEndpoint(1)).^2+curRotatedBallisticPredEndpoint(2).^2),'k+','linewidth',2,'markersize',10,'color',[0.65 0 0])
plot(371*ones(1,100), linspace(0,46.34,100),'k--','color',GOLD,'linewidth',1.5)
axis([0 endTime-startTime 0 inf])
xlabel('Time since go cue (ms)')
xticks([0 200 400 600])
% yticks([0 0.1 0.2 0.3 0.4])
ylabel('Distance from end target (mm)')
set(gca,'fontsize',10)

% 2D trajectory
subplot(4,2,[5 6])
plot(rotPos((firstInd:3:lastInd),1),rotPos((firstInd:3:lastInd),2),'k.-','markeredgecolor','k','markersize',10,'linewidth',0.5)
startTarget = 8.5*[cos(theta); sin(theta)]';
cursor = 3*[cos(theta); sin(theta)]';
hold on
plot(startTarget(:,1),startTarget(:,2),'k-','linewidth',0.5)
plot(target(:,1),target(:,2),'k-','linewidth',0.5)
plot(curRotatedBallisticPredEndpoint(1),curRotatedBallisticPredEndpoint(2),'k+','linewidth',2,'markersize',13,'color',[0.65 0 0])
plot(cursor(:,1)+curRotatedBallisticPredEndpoint(1),cursor(:,2)+curRotatedBallisticPredEndpoint(2),'k-','linewidth',0.5,'color',[0.65 0 0])
cursor = cursor+rotPos(lastInd,:);
plot(cursor(:,1),cursor(:,2),'k-','linewidth',0.5,'color',[0.65 0 0])
set(gca,'fontsize',10)
yticks([])
xticks([])
axis([-10 95 -20 20])


% Save it!
figname = 'Fig3A_BREPandHTSchematic';
saveas(gcf,['MATLABFigs\' figname])
saveas(gcf,['SVGs\' figname '.svg'])