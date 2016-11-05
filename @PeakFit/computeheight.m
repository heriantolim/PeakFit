function height=computeheight(peakShape,area,width)
%% Compute Peak Height
%  height=PeakFit.computeheight(peakShape,area,width) computes the height of a
%  peak function of the specified peakShape, given its areas and widths.
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 25/03/2013
% Last modified: 25/03/2013

assert(isintegerscalar(peakShape),...
	'PeakFit:computeheight:InvalidInput',...
	'Input to the peak shape must be an integer scalar.');
assert(isrealarray(area) && isrealarray(width),...
	'PeakFit:computeheight:InvalidInput',...
	'Input to the area and width must be arrays of real numbers.');
assert(all(size(area)==size(width)),...
	'PeakFit:computeheight:InvalidInput',...
	'The area and width arrays must have the same dimensions.');

switch peakShape
	case 1
		height=2/pi*area./width;
	case 2
		height=2/sqrt(pi/log(2))*area./width;
	otherwise
		error('PeakFit:computeheight:UnexpectedInput',...
			'The specified peak shape is not recognized.');
end

end