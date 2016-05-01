function [yModel,yPeak,yBaseline]=model(obj,varargin)
%% Fit Model
%  [yModel,yPeak,yBaseline]=obj.model() returns a set of fit models evaluated at
%  the data points.
%
%  [yModel,yPeak,yBaseline]=obj.model(x) evaluates the fit models at x instead.
%
% Outputs:
%  yModel: The reconstructed data points from the fit results.
%
%  yPeak: A cell containing the resconstructed data points of each peak.
%
%  yBaseline: The reconstructed baseline points.
%
% Requires package:
%  - Common_v1.0.0+
%
% Tested on:
%  - MATLAB R2013b
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 04/04/2013
% Last modified: 04/04/2013

yModel=[];
yPeak=[];
yBaseline=[];

if obj.NumPeaks==0
	return
end

%% Parse Inputs
numInputs=nargin-1;
switch numInputs
	case 0
		xModel=obj.XData;
	case 1
		assert(isrealvector(varargin{1}), ...
			'PeakFit:PeakFit:model:InvalidInput', ...
			['Input to the X data points for the model must be a vector of ',...
				'real numbers.']);
		xModel=varargin{1}(:)';
	otherwise
		error('PeakFit:PeakFit:model:UnexpectedInput',...
			'One or more inputs are not recognized.');
end
xModel=sort(xModel);

%% Initialization
numPeaks=obj.NumPeaks;
center=obj.Center;
height=obj.Height;
width=obj.Width;
baseline=obj.Baseline;

%% Construct Baseline
if obj.BaselinePolyOrder<0 || isempty(baseline)
	yBaseline=zeros(1,numel(xModel));
else
	yBaseline=polyval(baseline(1,:),xModel);
end

%% Construct Peak
yPeak=cell(1,numPeaks);
for i=1:numPeaks
	yPeak{i}=feval(@fnpeak,obj.PeakShape(i),xModel, ...
		center(1,i),height(1,i),width(1,i));
end

%% Construct Model
yModel=yBaseline;
for i=1:numPeaks
	yModel=yModel+yPeak{i};
end

end