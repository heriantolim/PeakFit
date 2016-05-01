function height=computepeakheight(peakShape,area,width)
%% Compute Peak Height
%  height=computepeakheight(peakShape,area,width) computes the height of a peak
%  function of the specified peakShape, given its area and width.
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 25/03/2013
% Last modified: 25/03/2013

assert(isstringscalar(peakShape), ...
	'PeakFit:computepeakheight:InvalidInput', ...
	'Input to the peak shape must be a string scalar.');
assert(isrealarray(area) && isrealarray(width), ...
	'PeakFit:computepeakheight:InvalidInput', ...
	'Input to the area and width must be arrays of real numbers.');
assert(all(size(area)==size(width)), ...
	'PeakFit:computepeakheight:InvalidInput', ...
	' The area and width vectors must have the same length.');

switch peakShape
	case 'Gaussian'
		height=2/sqrt(pi/log(2))*area./width;
	case 'Lorentzian'
		height=2/pi*area./width;
	otherwise
		error('PeakFit:computepeakheight:UnexpectedInput', ...
			'The specified peak shape is not recognized.');
end

end