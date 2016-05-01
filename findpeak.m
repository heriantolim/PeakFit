function [center,height]=findpeak(Data,varargin)
%% Find Peaks
%  [center,height]=findpeak(Data, ...) or
%  [center,height]=findpeak(xData,yData, ...) locates the maxima in the data
%  points. Data can be a m-by-2 or 2-by-n real matrix. xData and yData must be
%  real vectors of the same length.
%
%  A number of options can be specified to improve the accuracy of the peak
%  finding. The options, which are passed in a Name-Value pair argument format,
%  are:
%
%  MovingAvgWidth: To prevent identifying noise as peaks, the curve is smoothed 
%                  first by using simple moving average. By default, the width
%                  of the moving average is set to the ceil of the 1% of the
%                  number of coordinates. Set this options to override the
%                  default.
%
%  NumPeaks: Specify the number of peaks to be found. When more peaks are
%            found, the algorithm selects the peaks with the steepest
%            gradients. When NumPeaks=0, all peaks found are returned.
%            Defaults to 0.
%
%  Window: The algorithm only find peaks where the centers are within the
%          Window vector of length 2. When not specified, the algorithm
%          attempts to find all peaks in the given data points.
%
%  This function returns the centers and heights of the peaks found in a row
%  vector.
%
% Requires package:
%  - Common_v1.0.0+
%  - Math_v1.0.0+
%
% Tested on:
%  - MATLAB R2013b
%  - MATLAB R2015b
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 27/03/2013
% Last modified: 01/05/2016

%% Input Validation and Parsing
% Data points
numExtraInputs=nargin-1;
k=1;
if isempty(Data) || ~isrealmatrix(Data)
	error('PeakFit:findpeak:InvalidInput', ...
		'The first argument must be data points.');
elseif ~isvector(Data)
	[m,n]=size(Data);
	assert(m==2 || n==2, ...
		'PeakFit:findpeak:InvalidInput', ...
		'The matrix for the data points must be of size m-by-2 or 2-by-n.');
	if m~=2
		Data=Data';
		n=m;
	end
else
	assert(numExtraInputs>0, ...
		'PeakFit:findpeak:TooFewInputs', ...
		'The Y coordinates are required.');
	assert(isrealvector(varargin{1}), ...
		'PeakFit:findpeak:InvalidInput', ...
		['If the first argument is the X coordinates, ', ...
			'the second must be the Y coordinates.']);
	n=numel(Data);
	assert(numel(varargin{k})==n, ...
		'PeakFit:findpeak:InvalidInput', ...
		'The number of the X coordinates and Y coordinates must agree.');
	if iscolumn(Data)
		Data=Data';
	end
	if isrow(varargin{k})
		Data=[Data;varargin{k}];
	else
		Data=[Data;varargin{k}'];
	end
	k=k+1;
end
Data=sortcolumns(Data);

% Default options
movingAvgWidth=ceil(.01*n);
numPeaks=0;

% Name-Value pairs
while k<=numExtraInputs-1
	if ~isstringscalar(varargin{k})
		break;
	end
	switch varargin{k}
		case 'MovingAvgWidth'
			assert(isintegerscalar(varargin{k+1}) && varargin{k+1}>=0,...
				'PeakFit:findpeak:InvalidInput',...
				'Input to the MovingAvgWidth must be a positive integer.');
			movingAvgWidth=varargin{k+1};
		case 'NumPeaks'
			assert(isintegerscalar(varargin{k+1}) && varargin{k+1}>=0,...
				'PeakFit:findpeak:InvalidInput',...
				'Input to the NumPeaks must be a positive integer.');
			numPeaks=varargin{k+1};
		case 'Window'
			assert(isrealvector(varargin{k+1}) && numel(varargin{k+1})==2,...
				'PeakFit:findpeak:InvalidInput',...
				'Input to the Window must be a real vector of length two.');
			window=sort(varargin{k+1});
			Data=Data(:,Data(1,:)>=window(1) & Data(1,:)<=window(2));
			[~,n]=size(Data);
		otherwise
			break
	end
	k=k+2;
end
assert(k>numExtraInputs, ...
	'PeakFit:findpeak:UnexpectedInput', ...
	'One or more inputs are not recognized.');

%% Main
if numPeaks==1 % Special case
	[height,k]=max(Data(2,:));
	center=Data(1,k);
else
	numFoundPeaks=0;

	% Obtain the gradient of the curve
	Diff=diffquotient(Data(1,:),movingavg(Data(2,:),movingAvgWidth));

	% Find all maxima and give each a rank according to the average of their
	% adjacent absolute gradients
	for i=1:n-2
		if Diff(2,i)>0 && Diff(2,i+1)<0
			numFoundPeaks=numFoundPeaks+1;
			Rank(1,numFoundPeaks)=i+1; %#ok<AGROW>
			Rank(2,numFoundPeaks)=mean(abs(Diff(2,i:i+1))); %#ok<AGROW>
		end
	end
	
	if numFoundPeaks==0
		center=[];
		height=[];
	else
		if numPeaks>0 && numPeaks<numFoundPeaks
			% Sort the peaks in descending order based on their rank
			Rank=sortcolumns(Rank,-2);

			% Select only the best peaks, up to NumPeaks
			Rank=Rank(:,1:numPeaks);

			% Sort Rank back to its original order
			Rank=sortcolumns(Rank);

			numFoundPeaks=numPeaks;
		elseif numPeaks~=0
			warning('PeakFit:findpeak:FoundFewerPeaks', ...
				'Requested to find %d peaks, but only %d peaks are found.', ...
			numPeaks,numFoundPeaks);
		end
		
		center=zeros(1,numFoundPeaks);
		height=zeros(1,numFoundPeaks);
		for i=1:numFoundPeaks
			center(i)=Data(1,Rank(1,i));
			height(i)=Data(2,Rank(1,i));
		end
	end
end

end