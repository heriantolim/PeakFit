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
% Last modified: 25/10/2016

yModel=[];
yPeak=[];
yBaseline=[];

numPeaks=obj.NumPeaks;
center=obj.Center;
height=obj.Height;
width=obj.Width;

if numPeaks==0 || isempty(center) || isempty(height) || isempty(width)
	return
end

%% Parse Inputs
if nargin==1
	xModel=obj.XData;
elseif nargin>2
	error('PeakfFit:model:TooManyInput',...
		'At most one input argument is accepted.');
elseif isrealvector(varargin{1})
	xModel=varargin{1};
else
	error('PeakFit:model:InvalidInput',...
		'Input to the domains for the model must be a vector of real numbers.');
end

%% Construct Baseline
baseline=obj.Baseline;
if obj.BaselinePolyOrder<0 || isempty(baseline)
	yBaseline=zeros(size(xModel));
else
	yBaseline=polyval(baseline(1,:),xModel);
end

%% Construct Peak
yPeak=cell(1,numPeaks);
for i=1:numPeaks
	switch obj.PeakShape(i)
		case 1
			yPeak{i}=PeakFit.fnlorentzian(xModel,...
				center(1,i),height(1,i),width(1,i));
		case 2
			yPeak{i}=PeakFit.fngaussian(xModel,...
				center(1,i),height(1,i),width(1,i));
		otherwise
			error('PeakFit:model:UnexpectedInput',...
				'The specified peak shape is not recognized.');
	end
end

%% Construct Model
yModel=yBaseline;
for i=1:numPeaks
	yModel=yModel+yPeak{i};
end

end