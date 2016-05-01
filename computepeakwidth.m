function width=computepeakwidth(peakShape,area,height)
%% Compute Peak Width
%  width=computepeakwidth(peakShape,area,height) computes the width of a peak
%  function of the specified peakShape, given its area and height.
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 25/03/2013
% Last modified: 25/03/2013

assert(isstringscalar(peakShape), ...
	'PeakFit:computepeakwidth:InvalidInput', ...
	'Input to the peak shape must be a string scalar.');
assert(isrealarray(area) && isrealarray(height), ...
	'PeakFit:computepeakwidth:InvalidInput', ...
	'Input to the area and height must be arrays of real numbers.');
assert(all(size(area)==size(height)), ...
	'PeakFit:computepeakwidth:InvalidInput', ...
	'The area and height vectors must have the same length.');

switch peakShape
	case 'Gaussian'
		width=2/sqrt(pi/log(2))*area./height;
	case 'Lorentzian'
		width=2/pi*area./height;
	otherwise
		error('PeakFit:computepeakwidth:UnexpectedInput', ...
			'The specified peak shape is not recognized.');
end

end