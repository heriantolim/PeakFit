function [Diff,varargout]=diffquotient(Data,varargin)
%% (Forward) Difference Quotient
%  [Diff]=diffquotient(Data) or
%  [xDiff,yDiff]=diffquotient(xData,yData) calculates the difference quotient
%  (discrete derivative) to the first order from the data points. Data can be a
%  m-by-2 or 2-by-n real matrix; meanwhile xData and yData must be real vectors
%  of the same length.
%
%  [...]=diffquotient(...,order) calculates the difference quotient to the nth
%  order.
%
%  [...,yDiffError]=diffquotient(...,yDataError,order) also calculates the error
%  in the difference quotient, given that the yData has an error of yDataError.
%
%  [...,yDiffLow,yDiffUp]=diffquotient(...,yDataLow,yDataUp,order) also
%  calculates the error in the difference quotient, given that the yData has a
%  lower bound of yDataLow and an upper bound of yDataUp.
%
%  The outputs are sorted in ascending order based on the sorting of the xData;
%  and they preserve the shape of the input matrices.
%
% See also: diff, gradient.
%
% Requires package:
%  - Common_v1.0.0+
%
% Tested on:
%  - MATLAB R2013b
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 23/03/2013
% Last modified: 27/03/2013

%% Inputs Validation and Parsing
numExtraInputs=nargin-1;
k=1;

% Read the argument for data points
if isempty(Data) || ~isrealmatrix(Data)
	error('Math:diffquotient:InvalidInput', ...
		'The first argument must be data points.');
elseif ~isvector(Data)
	[m,n]=size(Data);
	assert(m==2 || n==2, ...
		'Math:diffquotient:InvalidInput', ...
		'The matrix for Data points must be of size m-by-2 or 2-by-n.');
	if m~=2
		Data=Data';
		n=m;
	end
	PARSE_DATA=1;
else    
	assert(numExtraInputs>0, ...
		'Math:diffquotient:TooFewInputs', ...
		'The Y coordinates are required.');
	assert(isrealvector(varargin{k}), ...
		'Math:diffquotient:InvalidInput', ...
		['If the first argument is the X coordinates, ', ...
			'the second must be the Y coordinates.']);
	[m,n]=size(Data);
	assert(numel(varargin{k})==m*n, ...
		'Math:diffquotient:InvalidInput', ...
		'The number of the X coordinates and Y coordinates must agree.');
	if iscolumn(Data)
		Data=Data';
	end
	if isrow(varargin{k})
		Data=[Data;varargin{k}];
	else
		Data=[Data;varargin{k}'];
	end
	PARSE_DATA=2;
	k=k+1;
end

% Read the argument for yDataError; OR yDataLow and yDataUp
PARSE_ERROR=0;
if k<=numExtraInputs
	if isrealvector(varargin{k}) && numel(varargin{k})==n
		PARSE_ERROR=1;
		if k<numExtraInputs
			if isrealvector(varargin{k+1}) && numel(varargin{k+1})==n
				PARSE_ERROR=2;
			end
		end
		switch PARSE_ERROR
			case 1
				if isrow(varargin{k})
					Data=[Data;varargin{k}];
				else
					Data=[Data;varargin{k}'];
				end
				k=k+1;
			case 2
				if isrow(varargin{k})
					Data=[Data;varargin{k}];
				else
					Data=[Data;varargin{k}'];
				end
				if isrow(varargin{k+1})
					Data=[Data;varargin{k+1}];
				else
					Data=[Data;varargin{k+1}'];
				end
				k=k+2;
		end
	end
end

% Read the argument for order
order=1;
if k<=numExtraInputs
	if isintegerscalar(varargin{k})
		order=varargin{k};
		k=k+1;
	end
end

assert(k>numExtraInputs, ...
	'Math:diffquotient:UnexpectedInput', ...
	'One or more inputs are not recognized.');

assert(n>order, ...
	'Math:diffquotient:InsufficientData', ...
	'The number of Data points must be greater than the order.');

%% Main
% Initialize output
numExtraOutputs=nargout-1;
k=1;

if PARSE_DATA+PARSE_ERROR>1
	varargout=cell(1,PARSE_DATA+PARSE_ERROR-1);
else
	varargout={};
end

% Sort
Data=sortcolumns(Data);

% Compute xDiff, yDiff
xDiff=mean([Data(1,1:(n-1));Data(1,2:n)]);
yDiff=diff(Data(2,:))./diff(Data(1,:));
switch PARSE_DATA
	case 1
		Diff=[xDiff;yDiff];
		if m~=2
			Diff=Diff';
		end
	case 2
		Diff=xDiff;
		if m~=1
			Diff=Diff';
		end
		if isrow(varargin{1})
			varargout{k}=yDiff;
		else
			varargout{k}=yDiff';
		end
		k=k+1;
end

% Compute yDiffError; OR yDiffLow and yDiffUp
switch PARSE_ERROR
	case 1
		yDiffError=mean([Data(3,1:(n-1));Data(3,2:n)]);
		if isrow(varargin{PARSE_DATA})
			varargout{k}=yDiffError;
		else
			varargout{k}=yDiffError';
		end
		k=k+1;
	case 2
		yDiffLow=(Data(3,2:n)-Data(4,1:(n-1)))./diff(Data(1,:));
		yDiffUp=(Data(4,2:n)-Data(3,1:(n-1)))./diff(Data(1,:));
		if isrow(varargin{PARSE_DATA})
			varargout{k}=yDiffLow;
		else
			varargout{k}=yDiffLow';
		end
		if isrow(varargin{PARSE_DATA+1})
			varargout{k+1}=yDiffUp;
		else
			varargout{k+1}=yDiffUp';
		end
		k=k+2;
end

% Prevent the outputs from being unassigned during call
assert(k>numExtraOutputs, ...
	'Math:diffquotient:TooManyOutputs', ...
	'The number of output assignments exceeds the limit.');

% Repeat process to compute the difference quotient of higher order  
for i=2:order
	[Diff,varargout{:}]=diffquotient(Diff,varargout{:});
end

end