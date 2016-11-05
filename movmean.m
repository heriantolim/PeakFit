function u=movmean(v,varargin)
%% Moving Mean
%  This function is provided for backwards compatibility. It does the same job
%  as the builtin function, movmean, which was introduced in MATLAB R2016a. It
%  can be safely deleted on newer MATLAB.
%
%  u=movmean(v) computes the simple moving average of the vector v using a
%  sample width=1.
%
%  u=movmean(v,w) computes the simple moving average of the vector v using a
%  sample width=w.
%
% Tested on:
%  - MATLAB R2013b
%  - MATLAB R2015b
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 22/03/2013
% Last modified: 20/10/2016

%% Input Validation and Parsing
assert(isvector(v) && isreal(v),...
	'Input to the sample vector must be a real vector.');
if isempty(varargin)
	w=1;
elseif isscalar(varargin{1}) && isreal(varargin{1}) && varargin{1}>=0
	w=varargin{1};
else
	error('Input to the sample width must be a positive real scalar.');
end

%% Simple Moving Average
u=v;
if w==0;
	return
end
n=numel(v);
w=ceil(w/2);
for i=2:w
	u(i)=mean(v(1:2*i-1));
	u(n-i+1)=mean(v(n-2*i+2:n));
end
n1=1+w;
n2=n-w;
for i=n1:n2
	u(i)=mean(v(i-w:i+w));
end

end