function obj=convertunit(obj,unitFrom,unitTo,varargin)
%% Convert Units
%  This method converts the values of the XData, Center, Baseline, and other
%  related properties in a PeakFit object according to a unit conversion.
%
%  Available units to be converted from or to are:
%    'eV'        : electron volt
%    'percm'     : cm^{-1} (wavenumber)
%    'Ramanshift': cm^{-1} (Raman shift)
%  When converting to or from Raman shift, an additional argument is required
%  that specifies the excitation wavelength in nano meter.
%
%  Examples:
%    obj=convertunit(obj,'nm','eV'); converts from nano meter to electron volt
%    obj=convertunit(obj,'eV','Ramanshift',532); converts from electron volt to
%    Raman shift at 532nm excitation.
%
% Requires package:
%  - Common_v1.0.0+
%  - PhysConst_v1.0.0+
%
% Tested on:
%  - MATLAB R2013b
%  - MATLAB R2015b
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 27/10/2016
% Last modified: 27/10/2016

%% Constants
DEF_UNITS={'eV','percm','Ramanshift'};

%% Input Validation and Parsing
assert(isstringscalar(unitFrom) && isstringscalar(unitTo),...
	'PeakFit:convertunit:InvalidInput',...
	'Input to the units name must be a string scalar.');
if strcmp(unitFrom,unitTo)
	return
end
assert(any(strcmp(unitFrom,DEF_UNITS)) && any(strcmp(unitTo,DEF_UNITS)),...
	'PeakFit:convertunit:UnexpectedCase',...
	'Conversion of the units from %s to %s has not been defined in the code.',...
	unitFrom,unitTo);

%% Conversion Coefficients
ME1=MException('PeakFit:convertunit:TooFewInput',...
	'Additional input arguments are required when converting from %s to %s.',...
	unitFrom,unitTo);
ME2=MException('PeakFit:convertunit:InvalidInput',...
	'The additional input arguments failed validation.');
switch unitFrom
	case 'eV'
		switch unitTo
			case 'percm'
				a=Constant.ElementaryCharge/Constant.Planck/Constant.LightSpeed/100;
				b=0;
			case 'Ramanshift'
				a=-Constant.ElementaryCharge/Constant.Planck/Constant.LightSpeed/100;
				if nargin<4
					throw(ME1);
				elseif ~isrealscalar(varargin{1}) && varargin{1}<=0
					throw(ME2);
				end
				b=1e7/varargin{1};
		end
	case 'percm'
		switch unitTo
			case 'eV'
				a=1e2*Constant.Planck*Constant.LightSpeed/Constant.ElementaryCharge;
				b=0;
			case 'Ramanshift'
				a=-1;
				if nargin<4
					throw(ME1);
				elseif ~isrealscalar(varargin{1}) && varargin{1}<=0
					throw(ME2);
				end
				b=1e7/varargin{1};
		end
	case 'Ramanshift'
		if nargin<4
			throw(ME1);
		elseif ~isrealscalar(varargin{1}) && varargin{1}<=0
			throw(ME2);
		end
		switch unitTo
			case 'eV'
				a=-1e2*Constant.Planck*Constant.LightSpeed/Constant.ElementaryCharge;
				b=-a*1e7/varargin{1};
			case 'percm'
				a=-1;
				b=1e7/varargin{1};
		end
end

%% Change the Object Properties
obj.XData=a*obj.XData+b;
obj.Window=a*obj.Window+b;

obj.CenterStart=a*obj.CenterStart+b;
if a<0
	obj.Center=a*obj.Center([1,3,2],:)+b;
	obj.CenterLow=a*obj.CenterUp+b;
	obj.CenterUp=a*obj.CenterLow+b;
else
	obj.Center=a*obj.Center+b;
	obj.CenterLow=a*obj.CenterLow+b;
	obj.CenterUp=a*obj.CenterUp+b;
end

if obj.BaselinePolyOrder>=0
	obj.Baseline=PeakFit.transformpolycoeff(obj.Baseline,a,b);
	B=[obj.BaselineStart;obj.BaselineLow;obj.BaselineUp];
	B=PeakFit.transformpolycoeff(B,a,b);
	obj.BaselineStart=B(1,:);
	obj.BaselineLow=B(2,:);
	obj.BaselineUp=B(3,:);
end

a=abs(a);
obj.Width=a*obj.Width;
obj.WidthStart=a*obj.WidthStart;
obj.WidthLow=a*obj.WidthLow;
obj.WidthUp=a*obj.WidthUp;
obj.AreaStart=a*obj.AreaStart;
obj.AreaLow=a*obj.AreaLow;
obj.AreaUp=a*obj.AreaUp;

end