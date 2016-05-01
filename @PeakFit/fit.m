function fit(obj)
%% Peak Fitting
%  This method performs a curve fitting to a spectral data with a linear
%  combination of peak functions using the MATLAB Curve Fitting toolbox. The fit
%  results which include the center, width, and height of the fitted peaks will
%  be stored in the object's properties. The confidence interval of the fit and
%  other fitting statistics will also be stored.
%
% Requires package:
%  - Math_v1.0.0+
%
% Tested on:
%  - MATLAB R2013b
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 15/03/2013
% Last modified: 02/05/2016

%% Constants
PEAK_SHAPE=1;
MOVING_AVG_WIDTH_FRAC=.01;
HEIGHT_START_FRAC=.8;
HEIGHT_LOW_FRAC=.05;

%% Initialization

% Put object properties to local memory
numPeaks=obj.NumPeaks;
baselinePolyOrder=obj.BaselinePolyOrder;
peakShape=obj.PeakShape;
areaStart=obj.AreaStart;
areaLow=obj.AreaLow;
areaUp=obj.AreaUp;
centerStart=obj.CenterStart;
centerLow=obj.CenterLow;
centerUp=obj.CenterUp;
heightStart=obj.HeightStart;
heightLow=obj.HeightLow;
heightUp=obj.HeightUp;
widthStart=obj.WidthStart;
widthLow=obj.WidthLow;
widthUp=obj.WidthUp;

% Check if data points are supplied
assert(~isempty(obj.XData) && ~isempty(obj.YData), ...
	'PeakFit:PeakFit:fit:MissingData', ...
	'The XData and YData are required to do a curve fitting.');
assert(numel(obj.XData)==numel(obj.YData), ...
	'PeakFit:PeakFit:fit:InconsistentNumPoints', ...
	'The number of X and Y data points must be the same.');

% Sort the points (x,y) in ascending order
Data=sortcolumns([obj.XData;obj.YData]);

% Update the object Data properties
obj.XData=Data(1,:);
obj.YData=Data(2,:);

% Trim data points
if ~isempty(obj.Window)
	Data=Data(:,Data(1,:)>=obj.Window(1) & Data(1,:)<=obj.Window(2));
else
	obj.Window=[Data(1,1),Data(1,end)];
end

% Initialize moving average width
if isempty(obj.MovingAvgWidth)
	obj.MovingAvgWidth=ceil(MOVING_AVG_WIDTH_FRAC*numel(Data)/2);
end

% Find min(y) and max(y)
yMax=max(Data(2,:));
if baselinePolyOrder<0
	yMin=0;
else
	yMin=min(Data(2,:));
end

% Perform linear mapping [x(1),x(end)] -> [0,1] and [min(y),max(y)] -> [0,1]
a=1/(Data(1,end)-Data(1,1));
b=-a*Data(1,1);
c=1/(yMax-yMin);
d=-c*yMin;
Data(1,:)=a*Data(1,:)+b;
Data(2,:)=c*Data(2,:)+d;
if ~isempty(areaStart)
	areaStart=a*c*areaStart;
end
if ~isempty(areaLow)
	areaLow=a*c*areaLow;
end
if ~isempty(areaUp)
	areaUp=a*c*areaUp;
end
if ~isempty(centerStart)
	centerStart=a*centerStart+b;
end
if ~isempty(centerLow)
	centerLow=a*centerLow+b;
end
if ~isempty(centerUp)
	centerUp=a*centerUp+b;
end
if ~isempty(heightStart)
	heightStart=c*heightStart;
end
if ~isempty(heightLow)
	heightLow=c*heightLow;
end
if ~isempty(heightUp)
	heightUp=c*heightUp;
end
if ~isempty(widthStart)
	widthStart=a*widthStart;
end
if ~isempty(widthLow)
	widthLow=a*widthLow;
end
if ~isempty(widthUp)
	widthUp=a*widthUp;
end

% Find the number of peaks
if numPeaks==0 && ...
		isempty(areaStart) && ...
		isempty(areaLow) && ...
		isempty(areaUp) && ...
		isempty(centerStart) && ...
		isempty(centerLow) && ...
		isempty(centerUp) && ...
		isempty(heightStart) && ...
		isempty(heightLow) && ...
		isempty(heightUp) && ...
		isempty(widthStart) && ...
		isempty(widthLow) && ...
		isempty(widthUp)
	[centerStart,heightStart]=findpeak(x,y,'MovingAvgWidth',obj.MovingAvgWidth);
	heightStart=heightStart*HEIGHT_START_FRAC;
end
numPeaks=max([numPeaks, ...
	numel(peakShape), ...
	numel(areaStart), ...
	numel(areaLow), ...
	numel(areaUp), ...
	numel(centerStart), ...
	numel(centerLow), ...
	numel(centerUp), ...
	numel(heightStart), ...
	numel(heightLow), ...
	numel(heightUp), ...
	numel(widthStart), ...
	numel(widthLow), ...
	numel(widthUp)]);

% More constants
WIDTH_START=.5/sqrt(numPeaks);
WIDTH_LOW=2*mean(diff(Data(1,:)));
WIDTH_UP=2*WIDTH_START;

% Fill the blank peak shape with defaults
if numel(peakShape)==1
	peakShape=peakShape*ones(1,numPeaks);
else
	peakShape=[peakShape,PEAK_SHAPE*ones(1,numPeaks-numel(peakShape))];
end

% Fill the blank constraints with NaN
areaStart=[areaStart,nan(1,numPeaks-numel(areaStart))];
areaLow=[areaLow,nan(1,numPeaks-numel(areaLow))];
areaUp=[areaUp,nan(1,numPeaks-numel(areaUp))];

centerStart=[centerStart,nan(1,numPeaks-numel(centerStart))];
centerLow=[centerLow,nan(1,numPeaks-numel(centerLow))];
centerUp=[centerUp,nan(1,numPeaks-numel(centerUp))];

heightStart=[heightStart,nan(1,numPeaks-numel(heightStart))];
heightLow=[heightLow,nan(1,numPeaks-numel(heightLow))];
heightUp=[heightUp,nan(1,numPeaks-numel(heightUp))];

widthStart=[widthStart,nan(1,numPeaks-numel(widthStart))];
widthLow=[widthLow,nan(1,numPeaks-numel(widthLow))];
widthUp=[widthUp,nan(1,numPeaks-numel(widthUp))];

% If centerStart is known, fill the blank heightStart
for i=1:numPeaks
	if isfinite(centerStart(i)) && ~isfinite(heightStart(i))
		heightStart(i)=interp1(Data(1,:),Data(2,:),centerStart(i)) ...
			*HEIGHT_START_FRAC;
	end
end

% If centerLow and centerUp exist, fill the blank centerStart and heightStart
for i=1:numPeaks
	if ~isfinite(centerStart(i)) && ...
			isfinite(centerLow(i)) && isfinite(centerUp(i))
		[center,height]=findpeak(Data(1,:),Data(2,:), ...
			'NumPeaks',1, ...
			'Window',[centerLow(i),centerUp(i)]);
		if isempty(center) || isempty(height)
			center=centerLow(i)+(centerUp(i)-centerLow(i))*rand();
			height=interp1(Data(1,:),Data(2,:),center);
		end
		centerStart(i)=center;
		if ~isfinite(heightStart(i))
			heightStart(i)=height*HEIGHT_START_FRAC;
		end
	end
end

% Mark any sequence of the remaining blanks in the centerStart
startIX=[];
endIX=[];
numSequences=1;
tf=false;
for i=1:numPeaks
	if ~isfinite(centerStart(i))
		if ~tf
			startIX(numSequences)=i; %#ok<AGROW>
		end
		tf=true;
	else
		if tf
			endIX(numSequences)=i-1; %#ok<AGROW>
			numSequences=numSequences+1;
		end
		tf=false;
	end
end
if ~isfinite(centerStart(numPeaks))
	endIX(numSequences)=numPeaks;
else
	numSequences=numSequences-1;
end

% Replace the blanks in CenterStart and HeightStart with default values
if numSequences>0
	searchWindow=zeros(2,numSequences);
	searchWindow(1,1)=Data(1,1)-obj.TolX;
	searchWindow(2,numSequences)=Data(1,end)+obj.TolX;
	for j=2:numSequences
		searchWindow(1,j)=centerStart(startIX(j-1))+obj.TolX;
	end
	for j=1:numSequences-1
		searchWindow(2,j)=centerStart(endIX(j+1))-obj.TolX;
	end
	for j=1:numSequences
		[center,height]=findpeak(Data(1,:),Data(2,:), ...
			'MovingAvgWidth',obj.MovingAvgWidth, ...
			'NumPeaks',endIX(j)-startIX(j)+1, ...
			'Window',searchWindow(1:2,j));
		centerStart(startIX(j):endIX(j))=center;
		for i=startIX(j):endIX(j)
			if ~isfinite(heightStart(i)) && ...
					~(isfinite(areaStart(i)) && isfinite(widthStart(i)))
				heightStart(i)=height(i-startIX(j)+1)*HEIGHT_START_FRAC;
			end
		end
	end
end

% Compute the width or height from the area
for i=1:numPeaks
	if isfinite(areaStart(i))
		if isfinite(heightStart(i)) && ~isfinite(widthStart(i))
			widthStart(i)=computepeakwidth(peakShape(i), ...
				areaStart(i),heightStart(i));
		elseif isfinite(widthStart(i)) && ~isfinite(heightStart(i))
			heightStart(i)=computepeakheight(peakShape(i), ...
				areaStart(i),widthStart(i));
		end
	end
	if isfinite(areaLow(i))
		if isfinite(heightLow(i)) && ~isfinite(widthUp(i))
			widthUp(i)=computepeakwidth(peakShape(i),areaLow(i),heightLow(i));
		elseif isfinite(widthLow(i)) && ~isfinite(heightUp(i))
			heightUp(i)=computepeakheight(peakShape(i),areaLow(i),widthLow(i));
		end
	end
	if isfinite(areaUp(i))
		if isfinite(heightUp(i)) && ~isfinite(widthLow(i))
			widthLow(i)=computepeakwidth(peakShape(i),areaUp(i),heightUp(i));
		elseif isfinite(widthUp(i)) && ~isfinite(heightLow(i))
			heightLow(i)=computepeakheight(peakShape(i),areaUp(i),widthUp(i));
		end
	end
end

% Replace the remaining blanks with defaults
for i=1:numPeaks
	if ~isfinite(widthStart(i))
		widthStart(i)=WIDTH_START;
	end
	if ~isfinite(widthLow(i))
		widthLow(i)=WIDTH_LOW;
	end
	if ~isfinite(widthUp(i))
		widthUp(i)=WIDTH_UP;
	end
	if ~isfinite(centerLow(i))
		centerLow(i)=centerStart(i)-widthStart(i);
	end
	if ~isfinite(centerUp(i))
		centerUp(i)=centerStart(i)+widthStart(i);
	end
	if ~isfinite(heightLow(i))
		heightLow(i)=heightStart(i)*HEIGHT_LOW_FRAC/HEIGHT_START_FRAC;
	end
	if ~isfinite(heightUp(i))
		heightUp(i)=heightStart(i)/HEIGHT_START_FRAC;
	end
end

% Compute the area from the height and width
for i=1:numPeaks
	if ~isfinite(areaStart(i))
		areaStart(i)=computepeakarea(peakShape(i),heightStart(i),widthStart(i));
	end
	if ~isfinite(areaLow(i))
		areaLow(i)=computepeakarea(peakShape(i),heightLow(i),widthLow(i));
	end
	if ~isfinite(areaUp(i))
		areaUp(i)=computepeakarea(peakShape(i),heightUp(i),widthUp(i));
	end
end

startPoint=[centerStart,heightStart,widthStart];
lowerBound=[centerLow,heightLow,widthLow];
upperBound=[centerUp,heightUp,widthUp];

if baselinePolyOrder>=0
	C=mean(Data(2,Data(1,:)<=Data(1,1)+.01*(Data(1,end)-Data(1,1))));
	baselineStart=[zeros(1,baselinePolyOrder),C];
	baselineLow=[-inf(1,baselinePolyOrder),C-.5];
	baselineUp=[inf(1,baselinePolyOrder),C+.5];
	if baselinePolyOrder>=1
		baselineLow(baselinePolyOrder)=-1;
		baselineUp(baselinePolyOrder)=1;
	end
	startPoint=[startPoint,baselineStart];
	lowerBound=[lowerBound,baselineLow];
	upperBound=[upperBound,baselineUp];    
end

%% Construct Fit Expression and Coefficients
expr=cell(1,numPeaks);
for i=1:numPeaks
	switch peakShape(i)
		case 1
			expr{i}=sprintf('h%1$d*(w%1$d/2)^2/((x-c%1$d).^2+(w%1$d/2)^2)',i);
		case 2
			expr{i}=sprintf('h%1$d*exp(-4*log(2)*(x-c%1$d).^2/w%1$d^2)',i);
		otherwise
			error('PeakFit:PeakFit:fit:UnexpectedCase', ...
				'Something went wrong. Please review the code.');
	end
end
expr=strjoin(expr,'+');
if baselinePolyOrder>=0
	for i=1:baselinePolyOrder-1
		expr=[expr,sprintf('+p%d*x.^%d',i,baselinePolyOrder-i+1)]; %#ok<AGROW>
	end
	if baselinePolyOrder>0
		expr=[expr,sprintf('+p%d*x',baselinePolyOrder)];
	end
	expr=[expr,sprintf('+p%d',baselinePolyOrder+1)];
end

coeff=cell(1,3*numPeaks);
for i=1:numPeaks
	coeff{i}=sprintf('c%d',i);
	coeff{numPeaks+i}=sprintf('h%d',i);
	coeff{2*numPeaks+i}=sprintf('w%d',i);
end
if baselinePolyOrder>=0
	coeff=[coeff,cell(1,baselinePolyOrder+1)];
	for i=1:baselinePolyOrder+1
		coeff{3*numPeaks+i}=sprintf('p%d',i);
	end
end

%% MATLAB Curve Fitting ToolBox
FitOption=fitoptions('Method',obj.Method, ...
	'Robust',obj.Robust, ...
	'StartPoint',startPoint, ...
	'Lower',lowerBound, ...
	'Upper',upperBound, ...
	'Algorithm',obj.Algorithm, ...
	'DiffMaxChange',obj.DiffMaxChange, ...
	'DiffMinChange',obj.DiffMinChange, ...
	'MaxFunEvals',obj.MaxFunEvals, ...
	'MaxIter',obj.MaxIters, ...
	'TolFun',obj.TolFun, ...
	'TolX',obj.TolX, ...
	'Display','off');
FitType=fittype(expr, ...
	'coefficients',coeff, ...
	'options',FitOption);
warning('off','all');
[FitObj,gof,FitOutput]=fit(Data(1,:)',Data(2,:)',FitType);
warning('on','all');

%% Extract Fit Results
fitResult=[coeffvalues(FitObj);confint(FitObj)];
center=fitResult(:,1:numPeaks);
height=fitResult(:,numPeaks+1:2*numPeaks);
width=fitResult(:,2*numPeaks+1:3*numPeaks);

% Reverse mapping [0,1] -> [x(1),x(end)] and [0,1] -> [min(y),max(y)]
center=(center-b)/a;
height=height/c;
width=width/a;
areaStart=areaStart/a/c;
areaLow=areaLow/a/c;
areaUp=areaUp/a/c;
centerStart=(centerStart-b)/a;
centerLow=(centerLow-b)/a;
centerUp=(centerUp-b)/a;
heightStart=heightStart/c;
heightLow=heightLow/c;
heightUp=heightUp/c;
widthStart=widthStart/a;
widthLow=widthLow/a;
widthUp=widthUp/a;

% Transform the baseline
if baselinePolyOrder>=0
	baseline=fitResult(:,3*numPeaks+1:3*numPeaks+baselinePolyOrder+1);
	if baselinePolyOrder>0
		A=zeros(baselinePolyOrder+1);
		for i=0:baselinePolyOrder
			for j=0:i
				A(baselinePolyOrder-i+1,baselinePolyOrder-i+j+1)= ...
					A(baselinePolyOrder-i+1,baselinePolyOrder-i+j+1) ...
					+nchoosek(i,j)*a^(i-j)*b^(j);
			end
		end
		B=zeros(size(baseline));
		for j=1:baselinePolyOrder+1
			for i=1:j
				B(1,j)=B(1,j)+baseline(1,i)*A(i,j);
				if A(i,j)<0
					B(2,j)=B(2,j)+baseline(3,j)*A(i,j);
					B(3,j)=B(3,j)+baseline(2,j)*A(i,j);
				else
					B(2,j)=B(2,j)+baseline(2,j)*A(i,j);
					B(3,j)=B(3,j)+baseline(3,j)*A(i,j);
				end
			end
		end
		baseline=B;
	end
	baseline(:,baselinePolyOrder+1)=baseline(:,baselinePolyOrder+1)-d;
	baseline=baseline/c;
end

% Update the object properties
obj.NumPeaks=numPeaks;
obj.PeakShape=peakShape;
obj.AreaStart=areaStart;
obj.AreaLow=areaLow;
obj.AreaUp=areaUp;
obj.CenterStart=centerStart;
obj.CenterLow=centerLow;
obj.CenterUp=centerUp;
obj.HeightStart=heightStart;
obj.HeightLow=heightLow;
obj.HeightUp=heightUp;
obj.WidthStart=widthStart;
obj.WidthLow=widthLow;
obj.WidthUp=widthUp;
obj.Center=center;
obj.Height=height;
obj.Width=width;
if baselinePolyOrder>=0
	obj.Baseline=baseline;
end
obj.Sse=gof.sse;
obj.R2=gof.rsquare;
obj.AdjustedR2=gof.adjrsquare;
obj.Std=gof.rmse;
obj.FirstOrderOptimality=FitOutput.firstorderopt;
obj.Residuals=FitOutput.residuals;
obj.Jacobian=FitOutput.Jacobian;
obj.ExitFlag=FitOutput.exitflag;
obj.NumObservations=FitOutput.numobs;
obj.NumCoeffs=FitOutput.numparam;
obj.NumFunEvals=FitOutput.funcCount;
obj.NumIters=FitOutput.iterations;

end