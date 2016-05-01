function value=get(obj,varargin)
%% Get
%  value=obj.get('PropertyName1','PropertyName2',...) retrieves the values of
%  the object's properties with the names given. If more than one property is
%  queried, the output value will be a 1-by-n cell containing the property
%  values, ordered in the sequence of the input property names.
%
% Requires package:
%  - Common_v1.0.0+
%
% Tested on:
%  - MATLAB R2013b
%
% See also: set.
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 17/03/2013
% Last modified: 02/05/2016

numInputs=nargin-1;
assert(numInputs>0, ...
	'PeakFit:PeakFit:get:ArgsRequired', ...
	'Please specify one or more property names to be queried.');
k=1;
value=cell(1,numInputs);
while k<=numInputs    
	if ~isstringscalar(varargin{k})
		break
	end
	metaInfo=obj.findprop(varargin{k});
	if isempty(metaInfo)
		break
	elseif strcmp(metaInfo.GetAccess,'public')
		value{k}=obj.(varargin{k});
		k=k+1;
	else
		error('PeakFit:PeakFit:get:RestrictedAccess', ...
			'Public access to get the %s is prohibited.', ...
			varargin{k});
	end
end
assert(k>numInputs, ...
	'PeakFit:PeakFit:get:UnexpectedInput', ...
	'One or more inputs are not recognized.');

end