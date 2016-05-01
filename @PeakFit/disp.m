function disp(obj)
%% Display
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 08/04/2013
% Last modified: 08/04/2013

numPeaks=obj.NumPeaks;
baselinePolyOrder=obj.BaselinePolyOrder;

fprintf('\nFit Options:\n');
fprintf('\t%20s: %d\n','MovingAvgWidth',obj.MovingAvgWidth);
fprintf('\t%20s: %s\n','Method',obj.Method);
fprintf('\t%20s: %s\n','Robust',obj.Robust);
fprintf('\t%20s: %s\n','Algorithm',obj.Algorithm);
fprintf('\t%20s: %e\n','DiffMaxChange',obj.DiffMaxChange);
fprintf('\t%20s: %e\n','DiffMinChange',obj.DiffMinChange);
fprintf('\t%20s: %d\n','MaxFunEvals',obj.MaxFunEvals);
fprintf('\t%20s: %d\n','MaxIters',obj.MaxIters);
fprintf('\t%20s: %e\n','TolFun',obj.TolFun);
fprintf('\t%20s: %e\n','TolX',obj.TolX);

header=sprintf(['   Lower Bound   Start Point   Upper Bound   ', ...
	' Lower CI     Convergence    Upper CI  \n']);
fprintf('\nArea:\n\t Peak');
fprintf(header);
lower=obj.AreaLow;
start=obj.AreaStart;
upper=obj.AreaUp;
x=obj.Area;
for i=1:numPeaks
	fprintf('\t%5d   %11.5e   %11.5e   %11.5e   %11.5e   %11.5e   %11.5e\n', ...
		i,lower(i),start(i),upper(i),x(2,i),x(1,i),x(3,i));
end

fprintf('\nCenter:\n\t Peak');
fprintf(header);
lower=obj.CenterLow;
start=obj.CenterStart;
upper=obj.CenterUp;
x=obj.Center;
for i=1:numPeaks
	fprintf('\t%5d   %11.5e   %11.5e   %11.5e   %11.5e   %11.5e   %11.5e\n', ...
		i,lower(i),start(i),upper(i),x(2,i),x(1,i),x(3,i));
end

fprintf('\nHeight:\n\t Peak');
fprintf(header);
lower=obj.HeightLow;
start=obj.HeightStart;
upper=obj.HeightUp;
x=obj.Height;
for i=1:numPeaks
	fprintf('\t%5d   %11.5e   %11.5e   %11.5e   %11.5e   %11.5e   %11.5e\n', ...
		i,lower(i),start(i),upper(i),x(2,i),x(1,i),x(3,i));
end

fprintf('\nWidth:\n\t Peak');
fprintf(header);
lower=obj.WidthLow;
start=obj.WidthStart;
upper=obj.WidthUp;
x=obj.Width;
for i=1:numPeaks
	fprintf('\t%5d   %11.5e   %11.5e   %11.5e   %11.5e   %11.5e   %11.5e\n', ...
		i,lower(i),start(i),upper(i),x(2,i),x(1,i),x(3,i));
end

fprintf('\nBaseline:\n\tCoeff');
fprintf(header);
lower=obj.BaselineLow;
start=obj.BaselineStart;
upper=obj.BaselineUp;
x=obj.Baseline;
for i=1:baselinePolyOrder+1
	fprintf('\t%5s   %11.5e   %11.5e   %11.5e   %11.5e   %11.5e   %11.5e\n', ...
		sprintf('p%d',i),lower(i),start(i),upper(i),x(2,i),x(1,i),x(3,i));
end

fprintf('\nGoodness of Fit:\n');
fprintf('\t%20s: %s\n','Sse',obj.Sse);
fprintf('\t%20s: %s\n','R2',obj.R2);
fprintf('\t%20s: %s\n','AdjustedR2',obj.AdjustedR2);
fprintf('\t%20s: %s\n','Std',obj.Std);
fprintf('\t%20s: %s\n','FirstOrderOptimality',obj.FirstOrderOptimality);

fprintf('\nIteration Report:\n');
fprintf('\t%20s: %d\n','NumIters',obj.NumIters);
fprintf('\t%20s: %d\n','NumFunEvals',obj.NumFunEvals);
fprintf('\t%20s: %d\n','ExitFlag',obj.ExitFlag);

end