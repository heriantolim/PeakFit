function area=computepeakarea(peakShape,height,width)
%% Compute Peak Area
%  width=computepeakwidth(peakShape,height,width) computes the width of a peak
%  function of the specified peakShape, given its height and width.
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 25/03/2013
% Last modified: 25/03/2013

assert(isstringscalar(peakShape), ...
    'PeakFit:computepeakarea:InvalidInput', ...
    'Input to the peak shape must be a string scalar.');
assert(isrealarray(height) && isrealarray(width), ...
    'PeakFit:computepeakarea:InvalidInput', ...
    'Input to the height and width must be arrays of real numbers.');
assert(all(size(height)==size(width)), ...
    'PeakFit:computepeakarea:InvalidInput', ...
    ' The height and width vectors must have the same length.');

switch peakShape
	case 'Gaussian'
		area=sqrt(pi/log(2))/2*height.*width;
	case 'Lorentzian'
		area=pi/2*height.*width;
	otherwise
		error('PeakFit:computepeakarea:UnexpectedInput', ...
			'The specified peak shape is not recognized.');
end

end