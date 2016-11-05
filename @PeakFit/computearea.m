function area=computearea(peakShape,height,width)
%% Compute Peak Area
%  area=PeakFit.computearea(peakShape,height,width) computes the area of a peak
%  function of the specified peakShape, given its heights and widths.
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 25/03/2013
% Last modified: 25/03/2013

assert(isintegerscalar(peakShape),...
	'PeakFit:computearea:InvalidInput',...
	'Input to the peak shape must be an integer scalar.');
assert(isrealarray(height) && isrealarray(width),...
	'PeakFit:computearea:InvalidInput',...
	'Input to the height and width must be arrays of real numbers.');
assert(all(size(height)==size(width)),...
	'PeakFit:computearea:InvalidInput',...
	'The height and width arrays must have the same dimensions.');

switch peakShape
	case 1
		area=pi/2*height.*width;
	case 2
		area=sqrt(pi/log(2))/2*height.*width;
	otherwise
		error('PeakFit:computearea:UnexpectedInput',...
			'The specified peak shape is not recognized.');
end

end