function width=computewidth(peakShape,area,height)
%% Compute Peak Width
%  width=PeakFit.computewidth(peakShape,area,height) computes the width of a
%  peak function of the specified peakShape, given its area and height.
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 25/03/2013
% Last modified: 25/03/2013

assert(isintegerscalar(peakShape),...
	'PeakFit:computewidth:InvalidInput',...
	'Input to the peak shape must be an integer scalar.');
assert(isrealarray(area) && isrealarray(height),...
	'PeakFit:computewidth:InvalidInput',...
	'Input to the area and height must be arrays of real numbers.');
assert(all(size(area)==size(height)),...
	'PeakFit:computewidth:InvalidInput',...
	'The area and height arrays must have the same dimensions.');

switch peakShape
	case 1
		width=2/pi*area./height;
	case 2
		width=2/sqrt(pi/log(2))*area./height;
	otherwise
		error('PeakFit:computewidth:UnexpectedInput',...
			'The specified peak shape is not recognized.');
end

end