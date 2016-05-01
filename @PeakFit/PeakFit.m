classdef PeakFit < handle
%% PeakFit Class
%  Fit a spectral curve with a linear combination of symmetric peak functions,
%  i.e. Gaussian, Lorentzian, etc. The fitting may include a baseline function
%  that represents a 'background' contribution to the spectra.
%
% Constructor:
%  obj=peakfit('PropertyName1','PropertyValue1',...),
%  obj=peakfit(Data,'PropertyName1','PropertyValue1',...), or
%  obj=peakfit(XData,YData,'PropertyName1','PropertyValue1',...) creates a
%  PeakFit object with data points specified in {XData, YData} and properties
%  set using the additional arguments. A m-by-2 or 2-by-n matrix Data can be
%  used to specify XData and YData. When the data points are given, the returned
%  object will have fit method executed in the construction.
%
% Set-able and Get-able Properties:
%
%  XData: The X data points of the curve.
%
%  YData: The Y data points of the curve.
%
%  Window: A vector of length two [a,b] that limits the fitting to only the
%          data points whose X coordinates lies within [a,b], where a<b.
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
%  ...Start: A vector of initial values for the ... coefficient. The
%            default values are determined heuristically. To make certain
%            elements default, use NaN. They will be then replaced with the
%            default values upon fitting.
%
%  ...Low: A vector of lower bounds on the ... coefficient to be fitted.
%          The default values are determined heuristically. To make certain
%          elements default, use -Inf. They will be then replaced with the
%          default values upon fitting.
%
%  ...Up: A vector of upper bounds on the ... coefficient to be fitted.
%         The default values are determined heuristically. To make certain
%         elements default, use Inf. They will be then replaced with the default
%         values upon fitting.
%
%  Area...: The start points, lower, or upper bounds for the area.
%
%  Center...: The start points, lower, or upper bounds for the center. By
%             default, the start point for the center is determined by finding
%             the maxima of the curve. When there are more maximas than the
%             specified NumPeaks, then the best maximas as many as NumPeaks are
%             selected based on the magnitude of the adjacent gradient. By
%             default, the lower and upper bounds of the center are the start
%             point minus and plus the width respectively.
%
%  Height...: The start points, lower, or upper bounds for the height. By
%             default, the upper bounds of the height is the maximum height of
%             the nearby maximas. The default for the start point is 0.8 of the
%             upper bounds, and the default for the lower point is 0.05 of the
%             upper bounds.
%
%  Width...: The start points, lower, or upper bounds for the width. By default,
%            the start point for the width is the range in XData divided by
%            2*sqrt(NumPeaks). The default for the lower bound is the average of
%            the difference in XData. The default for the upper bound is twice
%            the default for the start point.
%
%  Baseline...: The start points, lower, or upper bounds for the baseline.
%
%  IndependentVar: The independent variable for the fit expression. Defaults to
%                  'x'.
%
%  MovingAvgWidth: An integer that specifies the sample width of the moving
%                  average to be used for smoothing the curve in order to filter
%                  noise before finding the maximas. This parameter is used only
%                  when CenterStart is not given.
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
% Get-able Properties:
%  Method: The method used for the fitting, which is 'NonLinearLeastSquares'.
%
%  Peak: A struct containing the fit results for each peak.
%
%  Base: A struct containing the fit results for the baseline.
%
%  Area, Center, Height, Width, Baseline: A 3-by-NumPeaks matrix to store the
%        fit results for the area, ... , respectively. The first row is the
%        converged values; the second row is the 95% CI lower bounds; and the
%        third row is the 95% CI upper bounds.
%
%  Sse: Sum of Squares Error of the fit results.
%
%  R2: Coefficient of determination of the fit results.
%
%  AdjustedR2: Degree-of-freedom adjusted coefficient of determination of the
%              fit results.
%
%  Std: Standard Deviation of the fit results.
%
%  FirstOrderOptimality: Measure of the first-order optimality (absolute maximum
%                        of the gradient components).
%
%  Residuals: Vector of residuals.
%
%  Jacobian: Jacobian matrix.
%
%  ExitFlag: Describes the exit condition of the algorithm. Positive flags
%            indicate convergence, within tolerances. Zero flags indicate that
%            the maximum number of function evaluations or iterations was
%            exceeded. Negative flags indicate that the algorithm did not
%            converge to a solution.
%
%  NumObservations: Number of observations (response values).
%
%  NumCoeffs: Number of coeffcients to fit.
%
%  NumFunEvals: Number of function evaluations.
%
%  NumIters: Number of iterations.
%
% Public Methods:
%  disp: Display the fit options, start point, lower bounds, upper bounds, and
%        the fit results if any.
%
%  set: Set the values of the object properties.
%
%  get: Retrieve the values of the queried properties.
%
%  clone: Return an exact copy of the PeakFit object.
%
%  fit: Perform curve fitting. The results will be stored in the properties.
%
%  model: Return the reconstructed data points (a model) using the fit results.
%
% Requires package:
%  - Common_v1.0.0+
%  - Math_v1.0.0+
%
% Tested on:
%  - MATLAB R2013b
%  - MATLAB R2015b
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 25/03/2013
% Last modified: 01/05/2016

%% Properties
properties
	XData
	YData
	Window
	NumPeaks
	BaselinePolyOrder=0;
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
	MovingAvgWidth
	Robust='off';
	Algorithm='Trust-Region';
	DiffMaxChange=.1;
	DiffMinChange=1e-8;
	MaxFunEvals=50000;
	MaxIters=10000;
	TolFun=1e-6;
	TolX=1e-6;
end

properties %(SetAccess=private)
	Center
	Height
	Width
	Baseline
end

properties (Dependent=true,SetAccess=private)
	Area
	Peak
	Base
end

properties (Constant=true)
	Method='NonlinearLeastSquares';
end

properties (SetAccess=private)
	Sse
	R2
	AdjustedR2
	Std
	FirstOrderOptimality
	Residuals
	Jacobian
	ExitFlag
	NumObservations=0;
	NumCoeffs=0;
	NumFunEvals=0;
	NumIters=0;
end

properties (Constant=true,GetAccess=protected)
	PEAK_SHAPE_TABLE=table( ...
		[1;2], ...
		{'L';'G'}, ...
		{'Lorentzian';'Gaussian'}, ...
		'VariableNames',{'ID','Initial','Name'});
	ROBUST_LIST={'on','off','LAR','Bisquare'};
	ALGORITHM_LIST={'Levenberg-Marquardt','Gauss-Newton','Trust-Region'};
end

%% Methods
methods
	% Constructor
	function obj=PeakFit(varargin)
		obj.set(varargin{:});
		try
			obj.fit;
		catch exception
			if ~strcmp(exception.identifier,'PeakFit:PeakFit:fit:MissingData')
				throw(exception);
			end
		end
	end
	
	% Get Methods
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
				x(:,i)=computepeakarea(peakShape{i},height(:,i),width(:,i));
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
			warning('PeakFit:PeakFit:getPeak:MissingData', ...
				'Curve fitting has not been performed, call the fit method first.');
		else
			numPeaks=obj.NumPeaks;
			for i=1:numPeaks
				x(i).Area.Value=area(1,i);
				x(i).Area.CI=area(2:3,i)';
				x(i).Center.Value=center(1,i);
				x(i).Center.CI=center(2:3,i)';
				x(i).Height.Value=height(1,i);
				x(i).Height.CI=height(2:3,i)';
				x(i).Width.Value=width(1,i);
				x(i).Width.CI=width(2:3,i)';
			end
		end
	end
	
	function x=get.Base(obj)
		x=struct();
		baseline=obj.Baseline;
		if ~isempty(baseline)
			baselinePolyOrder=obj.BaselinePolyOrder;
			for i=1:baselinePolyOrder+1
				x.(sprintf('p%d',i)).Value=baseline(1,i);
				x.(sprintf('p%d',i)).Value=baseline(2:3,i)';
			end
		end
	end
	
	% Set Methods
	function set.XData(obj,x)
		if isempty(x)
			obj.XData=[];
		elseif isrealvector(x)
			obj.XData=x(:)';
		else
			error('PeakFit:PeakFit:setXData:InvalidInput', ...
				'Input to set the XData must be a vector of real numbers.');
		end
	end
	
	function set.YData(obj,x)
		if isempty(x)
			obj.YData=[];
		elseif isrealvector(x)
			obj.YData=x(:)';
		else
			error('PeakFit:PeakFit:setYData:InvalidInput', ...
				'Input to set the YData must be a vector of real numbers.');
		end
	end
	
	function set.Window(obj,x)
		if isempty(x)
			obj.Window=[];
		elseif isrealvector(x) && numel(x)==2
			obj.Window=x(:)';
		else
			error('PeakFit:PeakFit:setWindow:InvalidInput', ...
				'Input to set the Window must be a real vector of length two.');
		end
	end
	
	function set.NumPeaks(obj,x)
		if isintegerscalar(x) && x>=0
			obj.NumPeaks=x;
		else
			error('PeakFit:PeakFit:setNumPeaks:InvalidInput', ...
				'Input to set the NumPeaks must be a positive integer scalar.');
		end
	end
	
	function set.BaselinePolyOrder(obj,x)
		if isintegerscalar(x)
			obj.BaselinePolyOrder=x;
		else
			error('PeakFit:PeakFit:setBaselinePolyOrder:InvalidInput', ...
				'Input to set the BaselinePolyOrder must be an integer scalar.');
		end
	end
	
	function set.PeakShape(obj,x)
		if isempty(x)
			obj.PeakShape=[];
		elseif isintegervector(x)
			obj.PeakShape=x(:)';
		else
			if isstringscalar(x)
				x={x};
			elseif ~isstringvector(x)
				error('PeakFit:PeakFit:setPeakShape:InvalidInput',...
					'Input to PeakShape must be a string scalar or vector.');
			end
			ME=MException('PeakFit:PeakFit:setPeakShape:UnexpectedInput', ...
				'The specified peak shape has not been defined in the code.');
			N=numel(x);
			y=zeros(1,N);
			for n=1:N
				if numel(x(n))==1
					tf=strcmpi(x(n),obj.PEAK_SHAPE_TABLE.Initial);
				else
					tf=strcmpi(x(n),obj.PEAK_SHAPE_TABLE.Name);
				end
				if any(tf)
					y(n)=obj.PEAK_SHAPE_TABLE.ID(tf);
				else
					throw(ME);
				end
			end
			obj.PeakShape=y;
		end
	end
	
	function set.AreaStart(obj,x)
		if isempty(x)
			obj.AreaStart=[];
		elseif isrealvector(x)
			obj.AreaStart=x(:)';
		else
			error('PeakFit:PeakFit:setAreaStart:InvalidInput',...
				'Input to set the AreaStart must be a vector of real numbers.');
		end
	end
    
	function set.AreaLow(obj,x)
		if isempty(x)
			obj.AreaLow=[];
		elseif isrealvector(x)
			obj.AreaLow=x(:)';
		else
			error('PeakFit:PeakFit:setAreaLow:InvalidInput',...
				'Input to set the AreaLow must be a vector of real numbers.');
		end
	end
	
	function set.AreaUp(obj,x)
		if isempty(x)
			obj.AreaUp=[];
		elseif isrealvector(x)
			obj.AreaUp=x(:)';
		else
			error('PeakFit:PeakFit:setAreaUp:InvalidInput',...
				'Input to set the AreaUp must be a vector of real numbers.');
		end
	end
	
	function set.CenterStart(obj,x)
		if isempty(x)
			obj.CenterStart=[];
		elseif isrealvector(x)
			obj.CenterStart=x(:)';
		else
			error('PeakFit:PeakFit:setCenterStart:InvalidInput',...
				'Input to set the CenterStart must be a vector of real numbers.');
		end
	end
    
	function set.CenterLow(obj,x)
		if isempty(x)
			obj.CenterLow=[];
		elseif isrealvector(x)
			obj.CenterLow=x(:)';
		else
			error('PeakFit:PeakFit:setCenterLow:InvalidInput',...
				'Input to set the CenterLow must be a vector of real numbers.');
		end
	end
	
	function set.CenterUp(obj,x)
		if isempty(x)
			obj.CenterUp=[];
		elseif isrealvector(x)
			obj.CenterUp=x(:)';
		else
			error('PeakFit:PeakFit:setCenterUp:InvalidInput',...
				'Input to set the CenterUp must be a vector of real numbers.');
		end
	end
	
	function set.HeightStart(obj,x)
		if isempty(x)
			obj.HeightStart=[];
		elseif isrealvector(x)
			obj.HeightStart=x(:)';
		else
			error('PeakFit:PeakFit:setHeightStart:InvalidInput',...
				'Input to set the HeightStart must be a vector of real numbers.');
		end
	end
    
	function set.HeightLow(obj,x)
		if isempty(x)
			obj.HeightLow=[];
		elseif isrealvector(x)
			obj.HeightLow=x(:)';
		else
			error('PeakFit:PeakFit:setHeightLow:InvalidInput',...
				'Input to set the HeightLow must be a vector of real numbers.');
		end
	end
	
	function set.HeightUp(obj,x)
		if isempty(x)
			obj.HeightUp=[];
		elseif isrealvector(x)
			obj.HeightUp=x(:)';
		else
			error('PeakFit:PeakFit:setHeightUp:InvalidInput',...
				'Input to set the HeightUp must be a vector of real numbers.');
		end
	end
	
	function set.WidthStart(obj,x)
		if isempty(x)
			obj.WidthStart=[];
		elseif isrealvector(x)
			obj.WidthStart=x(:)';
		else
			error('PeakFit:PeakFit:setWidthStart:InvalidInput',...
				'Input to set the WidthStart must be a vector of real numbers.');
		end
	end
    
	function set.WidthLow(obj,x)
		if isempty(x)
			obj.WidthLow=[];
		elseif isrealvector(x)
			obj.WidthLow=x(:)';
		else
			error('PeakFit:PeakFit:setWidthLow:InvalidInput',...
				'Input to set the WidthLow must be a vector of real numbers.');
		end
	end
	
	function set.WidthUp(obj,x)
		if isempty(x)
			obj.WidthUp=[];
		elseif isrealvector(x)
			obj.WidthUp=x(:)';
		else
			error('PeakFit:PeakFit:setWidthUp:InvalidInput',...
				'Input to set the WidthUp must be a vector of real numbers.');
		end
	end
	
	function set.MovingAvgWidth(obj,x)
		if isempty(x)
			obj.MovingAvgWidth=[];
		elseif isintegerscalar(x)
			obj.MovingAvgWidth=x;
		else
			error('PeakFit:PeakFit:setMovingAvgWidth:InvalidInput',...
				'Input to set the MovingAvgWidth must be an integer scalar.');
		end
	end
	
	function set.Robust(obj,x)
		ME=MException('PeakFit:PeakFit:setRobust:InvalidInput',...
			['Input to set the Robust must be either: ', ...
				strjoin(obj.ROBUST_LIST,', '),'.']);
		if isstringscalar(x)
			tf=strcmpi(x,obj.ROBUST_LIST);
			if any(tf)
				obj.Robust=obj.ROBUST_LIST{tf};
			else
				throw(ME);
			end
		else
			throw(ME);
		end
	end
	
	function set.Algorithm(obj,x)
		ME=MException('PeakFit:PeakFit:setAlgorithm:InvalidInput',...
			['Input to set the Algorithm must be either: ', ...
				strjoin(obj.ALGORITHM_LIST,', '),'.']);
		if isstringscalar(x)
			tf=strcmpi(x,obj.ALGORITHM_LIST);
			if any(tf)
				obj.Algorithm=obj.ALGORITHM_LIST{tf};
			else
				throw(ME);
			end
		else
			throw(ME);
		end
	end
	
	function set.DiffMaxChange(obj,x)
		if isrealscalar(x) && x>0
			obj.DiffMaxChange=x;
		else
			error('PeakFit:PeakFit:setDiffMaxChange:InvalidInput',...
				'Input to set the DiffMaxChange must be a positive real scalar.');
		end
	end
	
	function set.DiffMinChange(obj,x)
		if isrealscalar(x) && x>0
			obj.DiffMinChange=x;
		else
			error('PeakFit:PeakFit:setDiffMinChange:InvalidInput',...
				'Input to set the DiffMinChange must be a positive real scalar.');
		end
	end
	
	function set.MaxFunEvals(obj,x)
		if isintegerscalar(x) && x>0
			obj.MaxFunEvals=x;
		else
			error('PeakFit:PeakFit:setMaxFunEvals:InvalidInput',...
				'Input to set the MaxFunEvals must be a positive integer scalar.');
		end
	end
	
	function set.MaxIters(obj,x)
		if isintegerscalar(x) && x>0
			obj.MaxIters=x;
		else
			error('PeakFit:PeakFit:setMaxIters:InvalidInput',...
				'Input to set the MaxIters must be a positive integer scalar.');
		end
	end
	
	function set.TolFun(obj,x)
		if isrealscalar(x) && x>0
			obj.TolFun=x;
		else
			error('PeakFit:PeakFit:setTolFun:InvalidInput',...
				'Input to set the TolFun must be a positive real scalar.');
		end
	end
	
	function set.TolX(obj,x)
		if isrealscalar(x) && x>0
			obj.TolX=x;
		else
			error('PeakFit:PeakFit:setTolX:InvalidInput',...
				'Input to set the TolX must be a positive real scalar.');
		end
	end
	
	disp(obj)
	set(obj,varargin)
	value=get(obj,varargin)
	fit(obj)
	[yModel,yPeak,yBaseline]=model(obj,varargin)
end

methods (Static=true)
	obj2=clone(obj1)
end

end