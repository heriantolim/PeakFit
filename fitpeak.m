function [center,height,width,PeakFitObj]=fitpeak(Data,varargin)
%% Fit Peaks
%  [...]=fitpeak(Data, ...) performs peak fitting to a set of data points given
%  in the first argument with the parameters given in the following arguments.
%
%  This function is a shortcut to create a PeakFit object and extract the
%  important computational outputs, namely: center, height, and width.
%  Confidence interval of these results are discarded, but they can be found in
%  the properties of the PeakFit object. The created object is returned in the
%  fourth output.
%
% Tested on:
%  - MATLAB R2013b
%  - MATLAB R2015b
%
% See also: PeakFit.
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 26/03/2013
% Last modified: 26/03/2013

PeakFitObj=PeakFit(Data,varargin{:});
center=PeakFitObj.Center;
height=PeakFitObj.Height;
width=PeakFitObj.Width;

center=center(1,:);
height=height(1,:);
width=width(1,:);

end