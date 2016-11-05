function disp(obj)
%% Display Object Properties
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 08/04/2013
% Last modified: 25/10/2016

fprintf('\nFit Options:\n');
fprintf('\t%18s:    %s\n','Method',obj.Method);
fprintf('\t%18s:    %s\n','Robust',obj.Robust);
fprintf('\t%18s:    %s\n','Algorithm',obj.Algorithm);
fprintf('\t%18s:    %.2e\n','DiffMaxChange',obj.DiffMaxChange);
fprintf('\t%18s:    %.2e\n','DiffMinChange',obj.DiffMinChange);
fprintf('\t%18s:    %d\n','MaxFunEvals',obj.MaxFunEvals);
fprintf('\t%18s:    %d\n','MaxIters',obj.MaxIters);
fprintf('\t%18s:    %.2e\n','TolFun',obj.TolFun);
fprintf('\t%18s:    %.2e\n','TolX',obj.TolX);

exitFlag=obj.ExitFlag;
if ~isempty(exitFlag)
	numPeaks=obj.NumPeaks;
	header=['   Lower Bound   Start Point   Upper Bound',...
	        '      Lower CI   Convergence      Upper CI     RelUncert\n'];
	format=['   %11.4e   %11.4e   %11.4e',...
	        '   %11.4e   %11.4e   %11.4e   %11.4e\n'];
	
	fprintf('\nArea:\n\t Peak');
	fprintf(header);
	lower=obj.AreaLow;
	start=obj.AreaStart;
	upper=obj.AreaUp;
	x=obj.Area;
	for i=1:numPeaks
		fprintf('\t%5d',i);
		fprintf(format,lower(i),start(i),upper(i),...
			x(2,i),x(1,i),x(3,i),abs((x(3,i)-x(2,i))/x(1,i))/2);
	end
	
	fprintf('\nCenter:\n\t Peak');
	fprintf(header);
	lower=obj.CenterLow;
	start=obj.CenterStart;
	upper=obj.CenterUp;
	x=obj.Center;
	for i=1:numPeaks
		fprintf('\t%5d',i);
		fprintf(format,lower(i),start(i),upper(i),...
			x(2,i),x(1,i),x(3,i),abs((x(3,i)-x(2,i))/x(1,i))/2);
	end
	
	fprintf('\nHeight:\n\t Peak');
	fprintf(header);
	lower=obj.HeightLow;
	start=obj.HeightStart;
	upper=obj.HeightUp;
	x=obj.Height;
	for i=1:numPeaks
		fprintf('\t%5d',i);
		fprintf(format,lower(i),start(i),upper(i),...
			x(2,i),x(1,i),x(3,i),abs((x(3,i)-x(2,i))/x(1,i))/2);
	end
	
	fprintf('\nWidth:\n\t Peak');
	fprintf(header);
	lower=obj.WidthLow;
	start=obj.WidthStart;
	upper=obj.WidthUp;
	x=obj.Width;
	for i=1:numPeaks
		fprintf('\t%5d',i);
		fprintf(format,lower(i),start(i),upper(i),...
			x(2,i),x(1,i),x(3,i),abs((x(3,i)-x(2,i))/x(1,i))/2);
	end
	
	Q=obj.BaselinePolyOrder+1;
	if Q>0
		fprintf('\nBaseline:\n\tCoeff');
		fprintf(header);
		lower=obj.BaselineLow;
		start=obj.BaselineStart;
		upper=obj.BaselineUp;
		x=obj.Baseline;
		for i=1:Q
			fprintf('\t%5s',sprintf('p%d',i));
			fprintf(format,lower(i),start(i),upper(i),...
				x(2,i),x(1,i),x(3,i),abs((x(3,i)-x(2,i))/x(1,i))/2);
		end
	end
	
	fprintf('\nGoodness of Fit:\n');
	fprintf('\t%18s:    %.4g\n','RelStDev',obj.RelStDev);
	fprintf('\t%18s:    %.4g\n','CoeffDeterm',obj.CoeffDeterm);
	fprintf('\t%18s:    %.4g\n','AdjCoeffDeterm',obj.AdjCoeffDeterm);
	
	fprintf('\nIteration Report:\n');
	fprintf('\t%18s:    %d\n','NumFunEvals',obj.NumFunEvals);
	fprintf('\t%18s:    %d\n','NumIters',obj.NumIters);
	if exitFlag>0
		exitStatus='converged';
	elseif exitFlag<0
		exitStatus='diverged';
	elseif obj.NumFunEvals>obj.MaxFunEvals
		exitStatus='MaxFunEvals was exceeded';
	elseif obj.NumIters>obj.MaxIters
		exitStatus='MaxIters was exceeded';
	else
		exitStatus='fail';
	end
	fprintf('\t%18s:    %d (%s)\n','ExitFlag',exitFlag,exitStatus);
end

end