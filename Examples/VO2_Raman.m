%% PeakFit Example #2: VO2 Raman
% 
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 29/10/2018
% Last modified: 29/10/2018

% Add the required packages using MatVerCon.
% addpackage('MatCommon','MatGraphics','PeakFit');

% Clear workspace variables.
clear;

% Load data. If fails, adjust the file path supplied to the argument.
S=load('VO2_Raman.mat');

% Perform Lorentzian peak fitting.
S.Fit=PeakFit(S.Data,'Window',[100,900],'PeakShape','Lorentzian',...
	'CenterLow',[521,145,198,224,258,304,310,336,387,394,430,503,588,615,665,827],...
	'CenterUp', [524,151,202,228,263,306,315,342,394,402,450,506,597,625,675,837],...
	'WidthUp',  [  4, 14, 10,  8, 10, 10, 10, 10, 18, 20, 40, 20, 30, 40, 20, 30],...
	'BaselinePolyOrder',1);

%% Plotting
% Settings.
Groot.usedefault();
Groot.usedefault('latex',8,.6);
RESOLUTION=300;
AXES_SIZE=[12,4];
TICK_LENGTH=.2;
CLIP_RANGE=[500,540];

% Plot data.
xData=S.Fit.XData;
yData=S.Fit.YData;
xLim=S.Fit.Window;
ix=(xData>=xLim(1)&xData<=xLim(2)) & (xData<CLIP_RANGE(1)|xData>CLIP_RANGE(2));
xModel=linspace(xLim(1),xLim(2),ceil(RESOLUTION/2.54*AXES_SIZE(1)));
[yModel,yPeak,yBaseline]=S.Fit.model(xModel);
yLim=[min(min(yData),min(yBaseline)),max(yData(ix))/.6];

% Figure.
fig=docfigure(AXES_SIZE);

% Axes.
pos=[0,0,AXES_SIZE];
ax=axes('Position',pos,'XLim',xLim,'YLim',yLim,'YTickLabel','','XDir','reverse');
xlabel('Raman shift (cm$^{-1}$)');
ylabel('Intensity (arb. unit)');
fixticklength(.2);

% Plots.
N=S.Fit.NumPeaks;
h=plot(xData,yData,'Color','b');
for j=1:N
	plot(xModel,yBaseline+yPeak{j},'Color','r','LineWidth',.3);
end
h(2)=plot(xModel,yModel,'Color','g','LineWidth',.3);

% Peak labels.
h=Label.peak(S.Fit.Center(1,:),'StringFormat','%.1f',...
	'FormatSpec',[2,ones(1,N-1)],'FontColor',{[.3,.3,.3],[1,.5,0]},...
	'LineWidth',0,'MinYPos',.1,'MinYDist',.01,'PlotLine',h,...
	'ClipRange',CLIP_RANGE);

% Reconfigure the axes.
extent=vertcat(h.Extent);
yLim(2)=max(extent(:,2)+extent(:,4));
yLim(2)=yLim(1)+diff(yLim)/(1-TICK_LENGTH/AXES_SIZE(2));
ax.YLim=yLim;

% Reconfigure the layout.
margin=ax.TightInset+.1;
ax.Position=pos+[margin(1),margin(2),0,0];
pos=pos+[0,0,margin(1)+margin(3),margin(2)+margin(4)];
set(fig,{'Position','PaperPosition','PaperSize'},{pos,pos,pos(3:4)});

% Saving.
print(fig,'VO2_Raman.png','-dpng',sprintf('-r%d',RESOLUTION));
close(fig);
