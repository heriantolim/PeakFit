function [xm,ym]=findmaxima(obj,varargin)
%% Find Maxima
%  [xm,ym]=obj.findmaxima() locates the maxima in the curve data stored in the
%  object and returns the coordinates of the maxima.
%
%  [xm,ym]=obj.findmaxima(window) locates the maxima within the specified x
%  window. window is a real vector of length two.
%
%  [xm,ym]=obj.findmaxima(n) locates exactly n maxima.
%
%  [xm,ym]=obj.findmaxima(window,n) locates n maxima within the window.
%
%  [xm,ym]=obj.findmaxima(x,y,...) use the x and y vector as the data points
%  instead of the data stored in the object.
%
%  The curve is smoothed first prior to finding the maxima with a simple moving
%  average/mean. The MovMeanWidth property is used as the sample width of the
%  moving mean.
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
% First created: 27/03/2013
% Last modified: 05/11/2016

%% Input Validation and Parsing
N=nargin-1;
if N==0
	x=obj.XData;
	y=obj.YData;
	n=0;
else
	k=N;
	if isintegerscalar(varargin{k}) && varargin{k}>=0
		n=varargin{k};
		k=k-1;
	else
		n=0;
	end
	if k>0 && isrealvector(varargin{k}) && numel(varargin{k})==2
		w=sort(varargin{k});
		k=k-1;
	else
		w=[];
	end
	if k>1 && isrealvector(varargin{k}) && isrealvector(varargin{k-1})...
			&& numel(varargin{k})==numel(varargin{k-1})
		y=varargin{k};
		x=varargin{k-1};
		k=k-2;
	else
		x=obj.XData;
		y=obj.YData;
	end
	if k>0
		error('PeakFit:findmaxima:UnexpectedInput',...
			'One or more inputs are not recognized.');
	end
end

% Trim data points
if ~isempty(w)
	ix=x>=w(1) & x<=w(2);
	x=x(ix);
	y=y(ix);
end

%% Finding the Maxima
N=numel(x);
if N==0
	% No data, return empty
	xm=[];
	ym=[];
elseif n==1 || N<4
	% The tallest point is always the maxima
	[ym,k]=max(y);
	xm=x(k);
else
	% The gradient of the smoothed curve
	w=obj.MovMeanWidth;
	if w<1
		w=ceil(w*N);
	end
	y1=diff(movmean(y,w))./diff(x);

	% Find all maxima and give each a rank according to the average of their
	% adjacent absolute gradients
	k=N-2;
	ix=zeros(1,k);
	rank=zeros(1,k);
	m=0;
	for i=1:k
		if y1(i)>0 && y1(i+1)<0
			m=m+1;
			ix(m)=i+1;
			rank(m)=mean(abs(y1(i:i+1)));
		end
	end
	ix=ix(1:m);
	rank=rank(1:m);
	
	if m==0
		% y is either monotonically increasing or decreasing
		[ym,k]=max(y);
		xm=x(k);
	else
		if n>0
			if m>n
				[~,k]=sort(rank,'descend');
				ix=ix(k);
				ix=sort(ix(1:n));
			elseif m<n
				warning('PeakFit:find:FoundFewerPeaks',...
					'Requested to find %d peaks, but only %d peaks are found.',...
					n,m);
			end
		end
		xm=x(ix);
		ym=y(ix);
	end
end

end
