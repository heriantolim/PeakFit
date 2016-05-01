function u=movingavg(v,varargin)
%% Moving Average
%  u=movingavg(v) computes the simple moving average of the vector v using a
%  sample width=1.
%
%  u=movingavg(v,w) computes the simple moving average of the vector v using a
%  sample width=w.
%
%  u=movingavg(v,method,...) computes the moving average of v according to the
%  given method. Input to the method must be a string scalar, e.g. 'simple',
%  'weighted', 'exponential'. The arguments after method are used to evaluate
%  the moving average.
%
% Requires package:
%  - Common_v1.0.0+
%
% Tested on:
%  - MATLAB R2013b
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 22/03/2013
% Last modified: 26/03/2013

%% Input Validation and Parsing
assert(isrealvector(v), ...
	'Math:movingavg:InvalidInput', ...
	'Input to the sample vector must be a real vector.');
numExtraInputs=nargin-1;
if numExtraInputs==0
	method='simple';
	width=1;
elseif isstringscalar(varargin{1})
	method=varargin{1};
	switch method
		case 'simple'
			if numExtraInputs==1
				width=1;
			elseif numExtraInputs==2
				assert(isintegerscalar(varargin{2}) && varargin{2}>=0, ...
					'Math:movingavg:InvalidInput', ...
					'Input to the sample width must be a positive integer.');
			else
				error('Math:movingavg:TooManyInputs', ...
					'Too many input arguments provided.');
			end
		case 'weighted'

		case 'exponential'

		otherwise
			error('Math:movingavg:UnexpectedInput', ...
				'The specified moving average method is not recognized.');        
	end
elseif isintegerscalar(varargin{1}) && varargin{1}>=0
	method='simple';
	width=varargin{1};
else
	error('Math:movingavg:UnexpectedInput', ...
		'One or more inputs are not recognized.');
end

switch method
	case 'simple'
		%% Simple Moving Average
		u=v;
		if width==0;
			return
		else
			n=numel(v);
			for i=2:width
				u(i)=mean(v(1:2*i-1));
				u(n-i+1)=mean(v(n-2*i+2:n));
			end
			for i=1+width:n-width
				u(i)=mean(v(i-width:i+width));
			end
		end
	% case ...

	otherwise
		error('Math:movingavg:IncompleteCode', ...
			'Code is incomplete. Please revise.');
end

end