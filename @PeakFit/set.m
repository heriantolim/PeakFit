function set(obj,varargin)
%% Set
%  obj.set('PropertyName1',PropertyValue1,...) sets the object's property named
%  'PropertyName1' with value PropertyValue1, and so on.
%
%  obj.set(Data,'PropertyName1',PropertyValue1,...), or
%  obj.set(XData,YData,'PropertyName1',PropertyValue1,...) sets the XData and
%  YData properties of the object. Data must be an m-by-2 or 2-by-n matrix.
%
% Requires package:
%  - Common_v1.0.0+
%
% Tested on:
%  - MATLAB R2013b
%
% See also: get.
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 17/03/2013
% Last modified: 02/05/2016

numInputs=nargin-1;
if numInputs==0
	return
end

k=1;
if isrealmatrix(varargin{1}) && ~isvector(varargin{1})
	[m,n]=size(varargin{1});
	assert(m==2 || n==2,...
		'PeakFit:PeakFit:set:InvalidInput', ...
		'Input to the data matrix must be of size m-by-2 or 2-by-n.');
	if m>2
		obj.XData=varargin{1}(:,1);
		obj.YData=varargin{1}(:,2);
	else
		obj.XData=varargin{1}(1,:);
		obj.YData=varargin{1}(2,:);
	end
	k=2;
elseif numInputs>1 && isrealvector(varargin{1}) && isrealvector(varargin{2})
	obj.XData=varargin{1};
	obj.YData=varargin{2};
	k=3;
end
while k<numInputs
	if ~isstringscalar(varargin{k})
		break
	end
	metaInfo=obj.findprop(varargin{k});
	if isempty(metaInfo)
		break
	elseif strcmp(metaInfo.SetAccess,'public')
		obj.(varargin{k})=varargin{k+1};
		k=k+2;
	else
		error('PeakFit:PeakFit:set:RestrictedAccess', ...
			'Public access to set the %s is prohibited.', ...
		varargin{k});
	end
end
assert(k>numInputs, ...
	'PeakFit:PeakFit:set:UnexpectedInput', ...
	'One or more inputs are not recognized.');

end