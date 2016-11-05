classdef PeakFit
%% PeakFit Class
%  Fit a spectral curve with a linear combination of symmetric peak functions,
%  e.g. Gaussian, Lorentzian, etc. A baseline polynomial function may be
%  included to represent a 'background' contribution to the spectra.
%
% Constructions:
%  There are two ways to properly construct a PeakFit object.
%
%  The first is by creating an empty PeakFit object and specifying the fitting
%  parameters via property assignments as follows:
%  obj=PeakFit();% creates an empty PeakFit
%  obj.XData=[...];% specify the X data points
%  obj.YData=[...];% specify the Y data points
%  obj.CenterStart=[...];% specify the start points of the peaks' center
%  obj.WidthStart=[...];% specify the start points for the peaks' width
%
%  To perform the peak fitting, call another PeakFit with the prepared obj as
%  the input argument.
%  obj=PeakFit(obj);
%
%  The fit results will be stored in the object properties, which can be
%  displayed into the console by calling:
%  disp(obj);
%
%  The second is by specifying the data points and fitting parameters together
%  in one call as follows:
%  obj=PeakFit(Data,...% specify the XData and YData
%     Window,[...],...% specify the fit window
%     CenterLow,[...],...% specify the lower bound for the peaks' center
%     CenterUp,[...],...% specify the upper bound for the peaks' center
%     WidthUp,[...],...% specify the upper bound for the peaks' width
%     BaselinePolyOrder,[.]);% specify the order for the polynomial baseline.
%  The data points can also be specified in the following way:
%  obj=PeakFit(XData,YData,...);
%
%  When the data points are supplied in the construction of a PeakFit object,
%  the fitting is automatically performed. The data points are the only
%  mandatory inputs for a PeakFit object.
%  
%  In most cases, specifying the following parameters:
%  Window, CenterLow, CenterUp, WidthUp, and BaselinePolyOrder,
%  as given in the last example, is sufficient to enable an accurate fitting.
%
% Public Properties:
%  Data: The data points of the curve. An alternative to XData and YData. Data
%        must be a two-column or two-row matrix.
%
%  XData: The X data points of the curve.
%
%  YData: The Y data points of the curve. The peak fitting may not properly if
%         some Y data points are negative.
%
%  Window: A vector of length two [a,b] that limits the fitting to only the
%          data points whose X coordinates lies within [a,b].
%
%  NumPeaks: The number of peaks wished to be fitted. When the fitting peak
%            shape, start point, lower, or upper bound are set with vectors of
%            length greater than NumPeaks, then NumPeaks will be incremented to
%            adjust to the maximum length of these vectors. When the maximum
%            length of these vectors is less, then these vectors will be
%            expanded and filled with the default values. When NumPeaks=0 and
%            all the start point, lower, and upper are not set, then the program
%            attempts to fit all peaks found in the curve.
%
%  PeakShape: A string vector that specifies the peak shape of each peak.
%             Available PeakShape in this version are: 'Lorentzian' (1) and
%             'Gaussian' (2). PeakShape may be set with an integer, the initial
%             of these names, e.g. 'L' or 'G', or the full name, e.g.
%             'Lorentzian'. When the length of PeakShape is less than NumPeaks,
%             the remaining peaks will have a default PeakShape, which is
%             'Lorentzian'. If PeakShape contains only one element, then the
%             default value is the same as that element.
%
%  ...Start: A vector of initial values for the ... coefficients, where the
%            blanks are one of the {Area,Center,Width,Height,Baseline}.
%            The default values are determined heuristically. To make certain
%            elements default, use NaN. They will be then replaced with the
%            default values upon fitting.
%
%  ...Low: A vector of lower bounds on the ... coefficients, where the
%          blanks are one of the {Area,Center,Width,Height,Baseline}.
%          The default values are determined heuristically. To make certain
%          elements default, use -Inf. They will be then replaced with the
%          default values upon fitting.
%
%  ...Up: A vector of upper bounds on the ... coefficients, where the
%         blanks are one of the {Area,Center,Width,Height,Baseline}.
%         The default values are determined heuristically. To make certain
%         elements default, use Inf. They will be then replaced with the default
%         values upon fitting.
%
%  Robust: Specifies the robust linear least-squares fitting method to be used.
%          Avaliable values are 'on', 'off', 'LAR', or 'Bisquare'. The default
%          is 'off'. 'LAR' specifies the least absolute residual method and
%          'Bisquare' specifies the bisquare weights method. 'on' is equivalent
%          to 'Bisquare', the default method.
%
%  Algorithm: The algorithm used for the fitting procedure. Available values are
%             'Lavenberg-Marquardt', 'Gauss-Newton', or 'Trust-Region'. The
%             default is 'Trust-Region'.
%
%  MovMeanWidth: The window width of the moving average used for smoothing the
%                curve in order to filter noise before finding the maximas. This
%                parameter is used only when CenterStart is not given. The value
%                of this can be set as a positive integer which specifies the
%                width in terms of the number of data points, OR a real scalar
%                between [0,1] which specifies the width in terms of a fraction
%                of the total number of data points.
%
%  DiffMaxChange: The maximum change in coefficients for finite difference
%                 gradients. The default is 0.1.
%
%  DiffMinChange: The minimum change in coefficients for finite difference
%                 gradients. The default is 10^-8.
%
%  MaxFunEvals: The maximum number of evaluations of the model allowed. The
%               default is 50000.
%
%  MaxIters: The maximum number of iterations allowed for the fit. The default
%            is 10000.
%
%  TolFun: The termination tolerance on the model value. The default is 10^-6.
%
%  TolX: The termination tolerance on the coefficient values. The default is
%        10^-6.
%
% Read-only Properties:
%  Method: The method used for the fitting, which is 'NonLinearLeastSquares'.
%
%  Peak: A struct containing the fit results for each peak.
%
%  Base: A struct containing the fit results for the baseline.
%
%  Area, Center, Height, Width, Baseline: A 3-by-NumPeaks matrix to store the
%        fit results for the area,... , respectively. The first row is the
%        converged values; the second row is the 95% CI lower bounds; and the
%        third row is the 95% CI upper bounds.
%
%  RelStDev: Relative Standard Deviation of the fit results.
%
%  CoeffDeterm: The coefficient of determination of the fit results.
%
%  AdjCoeffDeterm: The degree-of-freedom adjusted coefficient of determination
%                  of the fit results.
%
%  NumFunEvals: Number of function evaluations.
%
%  NumIters: Number of iterations.
%
%  ExitFlag: Describes the exit condition of the algorithm. Positive flags
%            indicate convergence, within tolerances. Zero flags indicate that
%            the maximum number of function evaluations or iterations was
%            exceeded. Negative flags indicate that the algorithm did not
%            converge to a solution.
%
% Public Methods:
%  disp: Display the options, results, error and performance of the fitting.
%
%  model: Return the reconstructed data points (model) using the fit results.
%
%  convertunit: Convert the units of the data points and the fit results.
%               Available units to be converted from or to are:
%                'eV'        : electron volt
%                'percm'     : cm^{-1} (wavenumber)
%                'Ramanshift': cm^{-1} (Raman shift)
%               When converting to or from Raman shift, an additional argument
%               is required that specifies the excitation wavelength in nano
%               meter.
%
% Static Methods:
%  fnlorentzian: The lorentzian function.
%
%  fngaussian: The gaussian function.
%
% Requires package:
%  - Common_v1.0.0+
%  - PhysConst_v1.0.0+ (for the convertunit method only)
%
% Tested on:
%  - MATLAB R2013b
%  - MATLAB R2015b
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 25/03/2013
% Last modified: 04/11/2016

%% Properties
% Data, start points, and constraints for fitting
properties (Dependent=true)
	Data
end
properties
	XData
	YData
	Window
	NumPeaks=0;
	PeakShape
	AreaStart
	AreaLow
	AreaUp
	CenterStart
	CenterLow
	CenterUp
	HeightStart
	HeightLow
	HeightUp
	WidthStart
	WidthLow
	WidthUp
	BaselinePolyOrder=0;
end

% The start point and constraints for the baseline are set automatically
properties (SetAccess=protected)
	BaselineStart
	BaselineLow
	BaselineUp
end

% Fitting options
properties (Constant=true)
	Method='NonlinearLeastSquares';
end
properties
	Robust='off';
	Algorithm='Trust-Region';
	MovMeanWidth=.02;
	DiffMaxChange=.1;
	DiffMinChange=1e-8;
	MaxFunEvals=1e5;
	MaxIters=1e3;
	TolFun=1e-6;
	TolX=1e-6;
end

% Fit results
properties (SetAccess=protected)
	Center
	Height
	Width
	Baseline
end
properties (Dependent=true)
	Area
	Peak
	Base
end

% Fit errors and performance
properties (SetAccess=protected)
	RelStDev=0;
	CoeffDeterm=0;
	AdjCoeffDeterm=0;
	NumFunEvals=0;
	NumIters=0;
	ExitFlag
end

% Class dictionary
properties (Constant=true,GetAccess=protected)
	DefaultPeakShape=1;
	MinNumPoints=10;
	PeakShapeTable=table(...
		[1;2],...
		{'L';'G'},...
		{'Lorentzian';'Gaussian'},...
		'VariableNames',{'ID','Initial','Name'});
	RobustList={'on','off','LAR','Bisquare'};
	AlgorithmList={'Levenberg-Marquardt','Trust-Region'};
end

%% Methods
methods
	% Constructor
	function obj=PeakFit(varargin)
		N=nargin;
		if N==0
			return
		end
		
		% Copy the PeakFit object given in the first argument if any
		k=1;
		if isa(varargin{k},class(obj))
			obj=varargin{k};
			k=k+1;
		end
		
		% Parse input to Data points
		if k<=N && isrealmatrix(varargin{k})
			if isvector(varargin{k})
				if k<N && isrealvector(varargin{k+1})
					obj.XData=varargin{k};
					obj.YData=varargin{k+1};
					k=k+2;
				end
			else
				obj.Data=varargin{k};
				k=k+1;
			end
		end
		
		% Parse inputs to the parameters
		P=properties(obj);
		while k<N
			if ~isstringscalar(varargin{k})
				break
			end
			ix=strcmpi(varargin{k},P);
			if any(ix)
				obj.(P{ix})=varargin{k+1};
				k=k+2;
			else
				break
			end
		end
		assert(k>N,...
			'PeakFit:UnexpectedInput',...
			'One or more inputs are not recognized.');
		
		% Perform peak fitting
		obj=fit(obj);
	end

	% Get Methods
	function x=get.Data(obj)
		if numel(obj.XData)==numel(obj.YData)
			x=[obj.XData;obj.YData];
		else
			x=[];
		end
	end
	
	function x=get.Area(obj)
		peakShape=obj.PeakShape;
		height=obj.Height;
		width=obj.Width;
		if isempty(peakShape) || isempty(height) || isempty(width)
			x=[];
		else
			numPeaks=numel(peakShape);
			x=zeros(3,numPeaks);
			for i=1:numPeaks
				x(:,i)=PeakFit.computearea(peakShape(i),height(:,i),width(:,i));
			end
		end
	end

	function x=get.Peak(obj)
		x=struct();
		area=obj.Area;
		center=obj.Center;
		height=obj.Height;
		width=obj.Width;
		if isempty(area) || isempty(center) || isempty(height) || isempty(width)
			return
		end
		numPeaks=obj.NumPeaks;
		for i=1:numPeaks
			x(i).Area.Value=area(1,i);
			x(i).Area.CI=area(2:3,i).';
			x(i).Center.Value=center(1,i);
			x(i).Center.CI=center(2:3,i).';
			x(i).Height.Value=height(1,i);
			x(i).Height.CI=height(2:3,i).';
			x(i).Width.Value=width(1,i);
			x(i).Width.CI=width(2:3,i).';
		end
	end

	function x=get.Base(obj)
		x=struct();
		baseline=obj.Baseline;
		if ~isempty(baseline)
			P=obj.BaselinePolyOrder;
			Q=P+1;
			for i=1:Q
				x.(sprintf('p%d',i)).Value=baseline(1,i);
				x.(sprintf('p%d',i)).CI=baseline(2:3,i).';
			end
		end
	end

	% Set Methods
	function obj=set.Data(obj,x)
		if isempty(x)
			obj.XData=[];
			obj.YData=[];
			return
		end
		ME=MException('PeakFit:InvalidInput',...
			'Input to set the Data must be a real matrix of size [n,2] or [2,n].');
		if isrealmatrix(x)
			[m,n]=size(x);
			if m==2
				obj.XData=x(1,:);
				obj.YData=x(2,:);
			elseif n==2
				obj.XData=x(:,1);
				obj.YData=x(:,2);
			else
				throw(ME);
			end
		else
			throw(ME);
		end
	end
	
	function obj=set.XData(obj,x)
		if isempty(x)
			obj.XData=[];
		elseif isrealvector(x) && all(isfinite(x))
			if numel(x)<obj.MinNumPoints
				error('PeakFit:setXData:InsufficientNumPoints',...
					'Input vector to the XData must have at least %d elements.',...
					obj.MinNumPoints);
			else
				obj.XData=x(:).';
			end
		else
			error('PeakFit:setXData:InvalidInput',...
				'Input to set the XData must be a vector of finite real numbers.');
		end
	end

	function obj=set.YData(obj,x)
		if isempty(x)
			obj.YData=[];
		elseif isrealvector(x) && all(isfinite(x))
			if numel(x)<obj.MinNumPoints
				error('PeakFit:setYData:InsufficientNumPoints',...
					'Input vector to the YData must have at least %d elements.',...
					obj.MinNumPoints);
			else
				obj.YData=x(:).';
			end
		else
			error('PeakFit:setYData:InvalidInput',...
				'Input to set the YData must be a vector of finite real numbers.');
		end
	end

	function obj=set.Window(obj,x)
		if isempty(x)
			obj.Window=[];
		elseif isrealvector(x) && numel(x)==2
			obj.Window=sort(x(:)).';
		else
			error('PeakFit:setWindow:InvalidInput',...
				'Input to set the Window must be a real vector of length two.');
		end
	end

	function obj=set.NumPeaks(obj,x)
		if isintegerscalar(x) && x>=0
			obj.NumPeaks=x;
		else
			error('PeakFit:setNumPeaks:InvalidInput',...
				'Input to set the NumPeaks must be a positive integer scalar.');
		end
	end

	function obj=set.PeakShape(obj,x)
		if isempty(x)
			obj.PeakShape=[];
		elseif isintegervector(x)
			obj.PeakShape=x(:).';
		else
			if isstringscalar(x)
				x={x};
			elseif ~isstringvector(x)
				error('PeakFit:setPeakShape:InvalidInput',...
					'Input to PeakShape must be a string scalar or vector.');
			end
			ME=MException('PeakFit:setPeakShape:UnexpectedInput',...
				'The specified peak shape has not been defined in this class.');
			N=numel(x);
			y=zeros(1,N);
			for n=1:N
				if numel(x{n})==1
					tf=strcmpi(x{n},obj.PeakShapeTable.Initial);
				else
					tf=strcmpi(x{n},obj.PeakShapeTable.Name);
				end
				if any(tf)
					y(n)=obj.PeakShapeTable.ID(tf);
				else
					throw(ME);
				end
			end
			obj.PeakShape=y;
		end
	end

	function obj=set.AreaStart(obj,x)
		if isempty(x)
			obj.AreaStart=[];
		elseif isrealvector(x)
			obj.AreaStart=x(:).';
		else
			error('PeakFit:setAreaStart:InvalidInput',...
				'Input to set the AreaStart must be a vector of real numbers.');
		end
	end

	function obj=set.AreaLow(obj,x)
		if isempty(x)
			obj.AreaLow=[];
		elseif isrealvector(x)
			obj.AreaLow=x(:).';
		else
			error('PeakFit:setAreaLow:InvalidInput',...
				'Input to set the AreaLow must be a vector of real numbers.');
		end
	end

	function obj=set.AreaUp(obj,x)
		if isempty(x)
			obj.AreaUp=[];
		elseif isrealvector(x)
			obj.AreaUp=x(:).';
		else
			error('PeakFit:setAreaUp:InvalidInput',...
				'Input to set the AreaUp must be a vector of real numbers.');
		end
	end

	function obj=set.CenterStart(obj,x)
		if isempty(x)
			obj.CenterStart=[];
		elseif isrealvector(x)
			obj.CenterStart=x(:).';
		else
			error('PeakFit:setCenterStart:InvalidInput',...
				'Input to set the CenterStart must be a vector of real numbers.');
		end
	end

	function obj=set.CenterLow(obj,x)
		if isempty(x)
			obj.CenterLow=[];
		elseif isrealvector(x)
			obj.CenterLow=x(:).';
		else
			error('PeakFit:setCenterLow:InvalidInput',...
				'Input to set the CenterLow must be a vector of real numbers.');
		end
	end

	function obj=set.CenterUp(obj,x)
		if isempty(x)
			obj.CenterUp=[];
		elseif isrealvector(x)
			obj.CenterUp=x(:).';
		else
			error('PeakFit:setCenterUp:InvalidInput',...
				'Input to set the CenterUp must be a vector of real numbers.');
		end
	end

	function obj=set.HeightStart(obj,x)
		if isempty(x)
			obj.HeightStart=[];
		elseif isrealvector(x)
			obj.HeightStart=x(:).';
		else
			error('PeakFit:setHeightStart:InvalidInput',...
				'Input to set the HeightStart must be a vector of real numbers.');
		end
	end

	function obj=set.HeightLow(obj,x)
		if isempty(x)
			obj.HeightLow=[];
		elseif isrealvector(x)
			obj.HeightLow=x(:).';
		else
			error('PeakFit:setHeightLow:InvalidInput',...
				'Input to set the HeightLow must be a vector of real numbers.');
		end
	end

	function obj=set.HeightUp(obj,x)
		if isempty(x)
			obj.HeightUp=[];
		elseif isrealvector(x)
			obj.HeightUp=x(:).';
		else
			error('PeakFit:setHeightUp:InvalidInput',...
				'Input to set the HeightUp must be a vector of real numbers.');
		end
	end

	function obj=set.WidthStart(obj,x)
		if isempty(x)
			obj.WidthStart=[];
		elseif isrealvector(x) && all(isnan(x) | x>=0)
			obj.WidthStart=x(:).';
		else
			error('PeakFit:setWidthStart:InvalidInput',...
				['Input to set the WidthStart must be a vector of positive ',...
					'real numbers.']);
		end
	end

	function obj=set.WidthLow(obj,x)
		if isempty(x)
			obj.WidthLow=[];
		elseif isrealvector(x) && all(isnan(x) | x>=0)
			obj.WidthLow=x(:).';
		else
			error('PeakFit:setWidthLow:InvalidInput',...
				['Input to set the WidthLow must be a vector of positive real ',...
					'numbers.']);
		end
	end

	function obj=set.WidthUp(obj,x)
		if isempty(x)
			obj.WidthUp=[];
		elseif isrealvector(x) && all(isnan(x) | x>=0)
			obj.WidthUp=x(:).';
		else
			error('PeakFit:setWidthUp:InvalidInput',...
				['Input to set the WidthUp must be a vector of positive real ',...
					'numbers.']);
		end
	end

	function obj=set.BaselinePolyOrder(obj,x)
		if isintegerscalar(x)
			obj.BaselinePolyOrder=x;
		else
			error('PeakFit:setBaselinePolyOrder:InvalidInput',...
				'Input to set the BaselinePolyOrder must be an integer scalar.');
		end
	end

	function obj=set.BaselineStart(obj,x)
		if isempty(x)
			obj.BaselineStart=[];
		elseif isrealvector(x)
			obj.BaselineStart=x(:).';
		else
			error('PeakFit:setBaselineStart:InvalidInput',...
				'Input to set the BaselineStart must be a vector of real numbers.');
		end
	end

	function obj=set.BaselineLow(obj,x)
		if isempty(x)
			obj.BaselineLow=[];
		elseif isrealvector(x)
			obj.BaselineLow=x(:).';
		else
			error('PeakFit:setBaselineLow:InvalidInput',...
				'Input to set the BaselineLow must be a vector of real numbers.');
		end
	end

	function obj=set.BaselineUp(obj,x)
		if isempty(x)
			obj.BaselineUp=[];
		elseif isrealvector(x)
			obj.BaselineUp=x(:).';
		else
			error('PeakFit:setBaselineUp:InvalidInput',...
				'Input to set the BaselineUp must be a vector of real numbers.');
		end
	end
	
	function obj=set.Robust(obj,x)
		ME=MException('PeakFit:setRobust:InvalidInput',...
			'Input to set the Robust must be either: %s.',...
			strjoin(obj.RobustList,', '));
		if isstringscalar(x)
			tf=strcmpi(x,obj.RobustList);
			if any(tf)
				obj.Robust=obj.RobustList{tf};
			else
				throw(ME);
			end
		else
			throw(ME);
		end
	end

	function obj=set.Algorithm(obj,x)
		ME=MException('PeakFit:setAlgorithm:InvalidInput',...
			'Input to set the Algorithm must be either: %s.',...
			strjoin(obj.AlgorithmList,', '));
		if isstringscalar(x)
			tf=strcmpi(x,obj.AlgorithmList);
			if any(tf)
				obj.Algorithm=obj.AlgorithmList{tf};
			else
				throw(ME);
			end
		else
			throw(ME);
		end
	end

	function obj=set.MovMeanWidth(obj,x)
		ME=MException('PeakFit:setMovMeanWidth:InvalidInput',...
			['Input to set the MovMeanWidth must be a positive integer scalar ',...
				'or a real scalar between [0,1].']);
		if isrealscalar(x) && x>=0
			if isintegerscalar(x) || x<=1
				obj.MovMeanWidth=x;
			else
				throw(ME);
			end
		else
			throw(ME);
		end
	end

	function obj=set.DiffMaxChange(obj,x)
		if isrealscalar(x) && x>0
			obj.DiffMaxChange=x;
		else
			error('PeakFit:setDiffMaxChange:InvalidInput',...
				'Input to set the DiffMaxChange must be a positive real scalar.');
		end
	end

	function obj=set.DiffMinChange(obj,x)
		if isrealscalar(x) && x>0
			obj.DiffMinChange=x;
		else
			error('PeakFit:setDiffMinChange:InvalidInput',...
				'Input to set the DiffMinChange must be a positive real scalar.');
		end
	end

	function obj=set.MaxFunEvals(obj,x)
		if isintegerscalar(x) && x>0
			obj.MaxFunEvals=x;
		else
			error('PeakFit:setMaxFunEvals:InvalidInput',...
				'Input to set the MaxFunEvals must be a positive integer scalar.');
		end
	end

	function obj=set.MaxIters(obj,x)
		if isintegerscalar(x) && x>0
			obj.MaxIters=x;
		else
			error('PeakFit:setMaxIters:InvalidInput',...
				'Input to set the MaxIters must be a positive integer scalar.');
		end
	end

	function obj=set.TolFun(obj,x)
		if isrealscalar(x) && x>0
			obj.TolFun=x;
		else
			error('PeakFit:setTolFun:InvalidInput',...
				'Input to set the TolFun must be a positive real scalar.');
		end
	end

	function obj=set.TolX(obj,x)
		if isrealscalar(x) && x>0
			obj.TolX=x;
		else
			error('PeakFit:setTolX:InvalidInput',...
				'Input to set the TolX must be a positive real scalar.');
		end
	end
	
	function obj=set.Center(obj,x)
		if ~isempty(x)
			if isrealmatrix(x) && size(x,1)==3
				obj.Center=x;
			else
				error('PeakFit:setCenter:InvalidInput',...
					'Input to set the Center must be a real matrix with 3 rows.');
			end
		end
	end
	
	function obj=set.Height(obj,x)
		if ~isempty(x)
			if isrealmatrix(x) && size(x,1)==3
				obj.Height=x;
			else
				error('PeakFit:setHeight:InvalidInput',...
					'Input to set the Height must be a real matrix with 3 rows.');
			end
		end
	end
	
	function obj=set.Width(obj,x)
		if ~isempty(x)
			if isrealmatrix(x) && size(x,1)==3
				obj.Width=x;
			else
				error('PeakFit:setWidth:InvalidInput',...
					'Input to set the Width must be a real matrix with 3 rows.');
			end
		end
	end
	
	function obj=set.Baseline(obj,x)
		if ~isempty(x)
			if isrealmatrix(x) && size(x,1)==3
				obj.Baseline=x;
			else
				error('PeakFit:setBaseline:InvalidInput',...
					'Input to set the Baseline must be a real matrix with 3 rows.');
			end
		end
	end
	
	function obj=set.RelStDev(obj,x)
		if isrealscalar(x) && x>=0
			obj.RelStDev=x;
		else
			error('PeakFit:setRelStDev:InvalidInput',...
				'Input to set the RelStDev must be a positive real scalar.');
		end
	end
	
	function obj=set.CoeffDeterm(obj,x)
		if isrealscalar(x)
			obj.CoeffDeterm=x;
		else
			error('PeakFit:setCoeffDeterm:InvalidInput',...
				'Input to set the CoeffDeterm must be a real scalar.');
		end
	end
	
	function obj=set.AdjCoeffDeterm(obj,x)
		if isrealscalar(x)
			obj.AdjCoeffDeterm=x;
		else
			error('PeakFit:setAdjCoeffDeterm:InvalidInput',...
				'Input to set the AdjCoeffDeterm must be a real scalar.');
		end
	end
	
	function obj=set.NumFunEvals(obj,x)
		if isintegerscalar(x) && x>=0
			obj.NumFunEvals=x;
		else
			error('PeakFit:setNumFunEvals:InvalidInput',...
				'Input to set the NumFunEvals must be a positive integer scalar.');
		end
	end
	
	function obj=set.NumIters(obj,x)
		if isintegerscalar(x) && x>=0
			obj.NumIters=x;
		else
			error('PeakFit:setNumIters:InvalidInput',...
				'Input to set the NumIters must be a positive integer scalar.');
		end
	end
	
	function obj=set.ExitFlag(obj,x)
		if isempty(x)
			obj.ExitFlag=[];
		elseif isintegerscalar(x)
			obj.ExitFlag=x;
		else
			error('PeakFit:setExitFlag:InvalidInput',...
				'Input to set the ExitFlag must be an integer scalar.');
		end
	end
	
	% display the object properties
	disp(obj)
	
	% construct a fit model from the fit results
	[yModel,yPeak,yBaseline]=model(obj,varargin)
	
	% convert units
	obj=convertunit(obj,unitFrom,unitTo,varargin)
end

methods (Access=protected)
	% fit peaks
	obj=fit(obj)
	
	% find maxima
	[xm,ym]=findmaxima(obj,varargin)
end

methods (Static=true)
	% lorentzian function
	y=fnlorentzian(x,c,h,w)
	
	% gaussian function
	y=fngaussian(x,c,h,w)
end

methods (Static=true, Access=protected)
	% compute the peak area
	area=computearea(peakShape,height,width)
	
	% compute the peak height
	height=computeheight(peakShape,area,width)
	
	% compute the peak width
	width=computewidth(peakShape,area,height)
	
	% transform polynomial coefficients
	p=transformpolycoeff(p,varargin)
end

end
