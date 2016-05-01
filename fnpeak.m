function y=fnpeak(peakShape,x,varargin)
%% Peak Function
%  y=fnpeak(peakShape,x, ...) evaluates a peak function corresponding to
%  peakShape at points x.
%
%  This function accepts extra arguments in Name-Value pairs, which will be
%  passed to the function.
%
% Requires package:
%  - Common_v1.0.0+
%
% Tested on:
%  - MATLAB R2013b
%  - MATLAB R2015b
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 04/04/2013
% Last modified: 01/05/2016

assert(isstringscalar(peakShape), ...
	'PeakFit:fnpeak:InvalidInput', ...
	'Input to the peak shape must be a string scalar.');
assert(isrealvector(x), ...
	'PeakFit:fnpeak:InvalidInput', ...
	'Input to x must be a vector of real numbers.');
if isintegerscalar(peakShape)
	switch peakShape
		case 1
			y=fnlorentzian(x,varargin{:});
		case 2
			y=fngaussian(x,varargin{:});
		otherwise
			error('PeakFit:fnpeak:UnexpectedInput', ...
				'The specified peak shape is not recognized.');
	end
elseif isstringscalar(peakShape)
	peakShape=lower(peakShape);
	switch peakShape
		case 'lorentzian'
			y=fnlorentzian(x,varargin{:});
		case 'gaussian'
			y=fngaussian(x,varargin{:});
		otherwise
			error('PeakFit:fnpeak:UnexpectedInput', ...
				'The specified peak shape is not recognized.');
	end
end

end