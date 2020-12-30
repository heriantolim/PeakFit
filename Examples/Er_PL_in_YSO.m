%% PeakFit Example #1: Er PL in YSO
%
% Copyright: Herianto Lim (https://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 28/10/2018
% Last modified: 28/10/2018

% Add the required packages using MatVerCon.
% addpackage('MatCommon','MatGraphics','PeakFit');

% Clear workspace variables.
clear;

% Load data. If fails, adjust the file path supplied to the argument.
S = load('Er_PL_in_YSO.mat');

% Perform Lorentzian peak fitting with default values.
% The algorithm will attempt to find all the peaks in the spectrum.
S.Fit = PeakFit(S.Data, 'PeakShape', 'Lorentzian');

%% Plotting
% Settings.
Groot.usedefault();
Groot.usedefault('latex', 8, .6);
RESOLUTION = 300;
AXES_SIZE = [12, 4];
TICK_LENGTH = .2;

% Plot data.
xData = S.Fit.XData;
yData = S.Fit.YData;
xLim = S.Fit.Window;
xModel = linspace(xLim(1), xLim(2), ceil(RESOLUTION / 2.54 * AXES_SIZE(1)));
[yModel, yPeak, yBaseline] = S.Fit.model(xModel);
yLim = [min(min(yData), min(yBaseline)), max(max(yData), max(yModel))];

% Figure.
fig = docfigure(AXES_SIZE);

% Axes.
pos = [0, 0, AXES_SIZE];
ax = axes('Position', pos, 'XLim', xLim, 'YLim', yLim, 'YTickLabel', '');
xlabel('Wavenumber (cm$^{-1}$)');
ylabel('Intensity (arb. unit)');
fixticklength(.2);

% Plots.
N = S.Fit.NumPeaks;
h = plot(xData, yData, 'Color', 'b');

for j = 1:N
    plot(xModel, yBaseline + yPeak{j}, 'Color', 'r', 'LineWidth', .3);
end

h(2) = plot(xModel, yModel, 'Color', 'g', 'LineWidth', .3);

% Peak labels.
h = Label.peak(S.Fit.Center(1, :), 'StringFormat', '%.1f', ...
    'FontColor', [.3, .3, .3], 'LineWidth', 0, ...
    'MinYPos', .1, 'MinYDist', .01, 'PlotLine', h);

% Reconfigure the axes.
extent = vertcat(h.Extent);
yLim(2) = max(extent(:, 2) + extent(:, 4));
yLim(2) = yLim(1) + diff(yLim) / (1 - TICK_LENGTH / AXES_SIZE(2));
ax.YLim = yLim;

% Reconfigure the layout.
margin = ax.TightInset + .1;
ax.Position = pos + [margin(1), margin(2), 0, 0];
pos = pos + [0, 0, margin(1) + margin(3), margin(2) + margin(4)];
set(fig, {'Position', 'PaperPosition', 'PaperSize'}, {pos, pos, pos(3:4)});

% Saving.
print(fig, 'Er_PL_in_YSO.png', '-dpng', sprintf('-r%d', RESOLUTION));
close(fig);
