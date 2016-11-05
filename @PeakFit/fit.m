function obj=fit(obj)
%% Peak Fitting
%  This method performs a curve fitting to a spectral data with a linear
%  combination of peak functions using the MATLAB Curve Fitting toolbox.
%
%  The fit results and statistics will be stored in the object properties after
%  a successful calling of this method.
%
%  The object properties related to the start points and constraints for the fit
%  parameters will also be initialized.
%
%  This algorithm is far from perfect. A lot of things need fixing to adapt to
%  different scenarios or be more efficient.
%
% Tested on:
%  - MATLAB R2013b
%  - MATLAB R2015b
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 15/03/2013
% Last modified: 06/11/2016

x=obj.XData;
y=obj.YData;
if isempty(x) || isempty(y)
	return
end

%% Defaults
HEIGHT_START=.85;
HEIGHT_UP=1.35;
WIDTH_START=.5;
WIDTH_UP=1;

%% Initialize the Fit Data
assert(numel(x)==numel(y),...
	'PeakFit:fit:InconsistentNumPoints',...
	'The number of X and Y data points must be equal.');

% Sort the points (x,y) in ascending order
[x,ix]=sort(x);
y=y(ix);

% Update the object data
obj.XData=x;
obj.YData=y;

% Trim data points
if ~isempty(obj.Window)
	ix=x>=obj.Window(1) & x<=obj.Window(2);
	assert(sum(ix)>=obj.MinNumPoints,...
		'PeakFit:fit:InsufficientNumPoints',...
		'The number of data points within the fit window must be at least %d.',...
		obj.MinNumPoints);
	x=x(ix);
	y=y(ix);
else
	obj.Window=[x(1),x(end)];
end

% Perform linear mapping from [x(1),x(end)]->[0,1] and [min(y),max(y)]->[0,1]
bp=obj.BaselinePolyOrder+1;
a=min(y);
b=max(y);
if bp<1 && a>0
	a=0;
end
c=1/(b-a);
d=-c*a;
a=1/(x(end)-x(1));
b=-a*x(1);
x=a*x+b;
y=c*y+d;

%% Initialize the Fit Parameters
% Find the number of peaks
n=[numel(obj.CenterStart),numel(obj.CenterLow),numel(obj.CenterUp),...
	numel(obj.HeightStart),numel(obj.HeightLow),numel(obj.HeightUp),...
	numel(obj.WidthStart),numel(obj.WidthLow),numel(obj.WidthUp),...
	numel(obj.AreaStart),numel(obj.AreaLow),numel(obj.AreaUp),...
	numel(obj.PeakShape),obj.NumPeaks];
np=max(n);

% Initialize the constraint and peak shape arrays
if np==0
	% Find all possible maxima
	[xm,ym]=obj.findmaxima(x,y);
	np=numel(xm);
	center=[xm;nan(2,np)];
	height=[HEIGHT_START*ym;nan(2,np)];
	width=nan(3,np);
	area=nan(3,np);
	ps=repmat(obj.DefaultPeakShape,1,np);
else
	% Fill the blank constraints with NaN
	n=np-n;
	center=[a*obj.CenterStart+b,nan(1,n(1));
		a*obj.CenterLow+b,nan(1,n(2));
		a*obj.CenterUp+b,nan(1,n(3))];
	height=[c*obj.HeightStart,nan(1,n(4));
		c*obj.HeightLow,nan(1,n(5));
		c*obj.HeightUp,nan(1,n(6))];
	width=[a*obj.WidthStart,nan(1,n(7));
		a*obj.WidthLow,nan(1,n(8));
		a*obj.WidthUp,nan(1,n(9))];
	area=[a*c*obj.AreaStart,nan(1,n(10));
		a*c*obj.AreaLow,nan(1,n(11));
		a*c*obj.AreaUp,nan(1,n(12))];
	if n(13)+1==np
		ps=repmat(obj.PeakShape,1,np);
	else
		ps=[obj.PeakShape,repmat(obj.DefaultPeakShape,1,n(13))];
	end
	
	% If lower bounds exceed upper bounds, swap them
	for j=1:np
		if center(2,j)>center(3,j)
			center([2,3],j)=center([3,2],j);
		end
		if height(2,j)>height(3,j)
			height([2,3],j)=height([3,2],j);
		end
		if width(2,j)>width(3,j)
			width([2,3],j)=width([3,2],j);
		end
		if area(2,j)>area(3,j)
			area([2,3],j)=area([3,2],j);
		end
	end
end

% If center starts are known, fill the blank height starts
for j=1:np
	if isfinite(center(1,j)) && ~isfinite(height(1,j))
		height(1,j)=HEIGHT_START*interp1(x,y,center(1,j));
	end
end

% If the lower and/or upper bounds exist, fill the blank start points
for j=1:np
	if all(isfinite(center(:,j))==[0;1;1])
		[center(1,j),ym]=obj.findmaxima(x,y,center(2:3,j),1);
		if ~isfinite(height(1,j))
			height(1,j)=HEIGHT_START*ym;
		end
	end
	if ~isfinite(height(1,j))
		if isfinite(height(3,j))
			if isfinite(height(2,j))
				height(1,j)=height(2,j)+HEIGHT_START*(height(3,j)-height(2,j));
			else
				height(1,j)=HEIGHT_START*height(3,j);
			end
		end
	end
	if ~isfinite(width(1,j))
		if isfinite(width(3,j))
			if isfinite(width(2,j))
				width(1,j)=width(2,j)+WIDTH_START*(width(3,j)-width(2,j));
			else
				width(1,j)=WIDTH_START*width(3,j);
			end
		end
	end
	if ~isfinite(area(1,j))
		if isfinite(area(3,j))
			if isfinite(area(2,j))
				area(1,j)=area(2,j)+HEIGHT_START*WIDTH_START*(area(3,j)-area(2,j));
			else
				area(1,j)=HEIGHT_START*WIDTH_START*area(3,j);
			end
		end
	end
end

% Fill the remaining blank center starts.
% This block needs revision because center starts may not be initially ordered.
ix=~isfinite(center(1,:));
if any(ix)
	ix=find(ix);
	w=diff(ix);
	s=ix([2,w]>1);
	t=ix([w,2]>1);
	n=numel(s);
	w=zeros(n,2);
	if s(1)>1
		w(1,1)=center(1,s(1)-1)+obj.TolX;
	else
		w(1,1)=x(1)-obj.TolX;
	end
	if t(n)<np
		w(n,2)=center(1,t(n)+1)-obj.TolX;
	else
		w(n,2)=x(end)+obj.TolX;
	end
	for i=2:n
		w(i,1)=center(1,s(i)-1)+obj.TolX;
		w(i-1,2)=center(1,t(i-1)+1)-obj.TolX;
	end
	for i=1:n
		[xm,ym]=obj.findmaxima(x,y,w(i,:),t(i)-s(i)+1);
		center(1,s(i):t(i))=xm;
		for j=s(i):t(i)
			if ~isfinite(height(1,j))...
					&& ~(isfinite(area(1,j)) && isfinite(width(1,j)))
				height(1,j)=HEIGHT_START*ym(j-s(i)+1);
			end
		end
	end
end

% Compute the width or height constraints from the area constraints
n=[1,3,2];
for i=1:3
	for j=1:np
		if isfinite(area(i,j))
			if ~isfinite(width(i,j)) && isfinite(height(n(i),j))
				width(i,j)=PeakFit.computewidth(ps(j),area(i,j),height(n(i),j));
			elseif ~isfinite(height(i,j)) && isfinite(width(n(i),j))
				if area(i,j)<0
					height(i,j)=PeakFit.computeheight(ps(j),area(n(i),j),width(i,j));
				else
					height(i,j)=PeakFit.computeheight(ps(j),area(i,j),width(n(i),j));
				end
			end
		end
	end
end

% Replace the blank width constraints with defaults
width(3,~isfinite(width(3,:)))=WIDTH_UP/sqrt(np);
width(2,~isfinite(width(2,:)))=2*mean(diff(x));
ix=~isfinite(width(1,:));
width(1,ix)=max(width(2,ix),WIDTH_START*width(3,ix));

% Replace the blank center constraints with defaults
ix=~isfinite(center(1,:));
center(1,ix)=rand(1,sum(ix));
ix=~isfinite(center(2,:));
center(2,ix)=center(1,ix)-width(1,ix);
ix=~isfinite(center(3,:));
center(3,ix)=center(1,ix)+width(1,ix);

% Replace the blank height constraints with defaults
ix=~isfinite(height(1,:));
height(1,ix)=rand(1,sum(ix));
w=2*obj.TolFun;
for j=1:np
	if height(1,j)<0
		if ~isfinite(height(2,j))
			height(2,j)=-HEIGHT_UP*min(y(x>=center(2,j) & x<=center(3,j)));
		end
		if ~isfinite(height(3,j))
			height(3,j)=-w;
		end
	else
		if ~isfinite(height(2,j))
			height(2,j)=w;
		end
		if ~isfinite(height(3,j))
			height(3,j)=HEIGHT_UP*max(y(x>=center(2,j) & x<=center(3,j)));
		end
	end
end

% Compute the area from the height and width
for i=1:3
	for j=1:np
		if ~isfinite(area(i,j))
			area(i,j)=PeakFit.computearea(ps(j),height(i,j),width(i,j));
		end
	end
end

% Set the baseline constraints
if bp>0
	% Set the y-axis intersect as the start point for the baseline
	baseline=zeros(1,bp);
	w=obj.MovMeanWidth;
	if w>=1
		w=w/numel(x);
	end
	baseline(bp)=mean(y(x<=w/2));
	
	% Use a fibonacci sequence to define the lower and upper bounds
	n=bp+1;
	w=ones(1,n);
	for i=3:n
		w(i)=w(i-1)+w(i-2);
	end
	w=flip(w);
	w(n)=[];
	baseline=[baseline;baseline-w;baseline+w];
else
	baseline=zeros(3,0);
end

%% Construct the Fit Coefficients and Expression
coeff=cell(1,3*np);
for i=1:np
	coeff{i}=sprintf('c%d',i);
	coeff{np+i}=sprintf('h%d',i);
	coeff{2*np+i}=sprintf('w%d',i);
end
if bp>0
	coeff=[coeff,cell(1,bp)];
	for i=1:bp
		coeff{3*np+i}=sprintf('p%d',i);
	end
end

expr=cell(1,np);
for i=1:np
	switch ps(i)
		case 1
			expr{i}=sprintf('h%1$d*(w%1$d/2)^2/((x-c%1$d).^2+(w%1$d/2)^2)',i);
		case 2
			expr{i}=sprintf('h%1$d*exp(-4*log(2)*(x-c%1$d).^2/w%1$d^2)',i);
		otherwise
			error('PeakFit:fit:UnexpectedCase',...
				'Something went wrong. Please review the code.');
	end
end
expr=strjoin(expr,'+');
if bp>0
	n=bp-2;
	for i=1:n
		expr=[expr,sprintf('+p%d*x.^%d',i,bp-i)]; %#ok<AGROW>
	end
	if bp>1
		expr=[expr,sprintf('+p%d*x',bp-1)];
	end
	expr=[expr,sprintf('+p%d',bp)];
end

%% MATLAB Curve Fitting ToolBox
FitOption=fitoptions('Method',obj.Method,...
	'Robust',obj.Robust,...
	'StartPoint',[center(1,:),height(1,:),width(1,:),baseline(1,:)],...
	'Lower',[center(2,:),height(2,:),width(2,:),baseline(2,:)],...
	'Upper',[center(3,:),height(3,:),width(3,:),baseline(3,:)],...
	'Algorithm',obj.Algorithm,...
	'DiffMaxChange',obj.DiffMaxChange,...
	'DiffMinChange',obj.DiffMinChange,...
	'MaxFunEvals',obj.MaxFunEvals,...
	'MaxIter',obj.MaxIters,...
	'TolFun',obj.TolFun,...
	'TolX',obj.TolX,...
	'Display','off');
FitType=fittype(expr,...
	'coefficients',coeff,...
	'options',FitOption);
warning('off','all');
[FitObj,gof,FitOutput]=fit(x',y',FitType);
warning('on','all');
fitResult=[coeffvalues(FitObj);confint(FitObj)];

%% Update Object Properties
% Update the shape and number of peaks
obj.NumPeaks=np;
obj.PeakShape=ps;

% Reverse mapping coefficients
a=1/a;
b=-a*b;
c=1/c;
d=-c*d;

% Update the center, height, and width, after a reverse mapping
center=a*center+b;
height=c*height;
width=a*width;
area=a*c*area;
obj.Center=a*fitResult(:,1:np)+b;
obj.Height=c*fitResult(:,np+1:2*np);
obj.Width=a*fitResult(:,2*np+1:3*np);
obj.CenterStart=center(1,:);
obj.CenterLow=center(2,:);
obj.CenterUp=center(3,:);
obj.HeightStart=height(1,:);
obj.HeightLow=height(2,:);
obj.HeightUp=height(3,:);
obj.WidthStart=width(1,:);
obj.WidthLow=width(2,:);
obj.WidthUp=width(3,:);
obj.AreaStart=area(1,:);
obj.AreaLow=area(2,:);
obj.AreaUp=area(3,:);

% Update the baseline, after a reverse mapping
if bp>0
	baseline=PeakFit.transformpolycoeff(baseline,a,b,c,d);
	obj.Baseline=PeakFit.transformpolycoeff(fitResult(:,3*np+1:3*np+bp),a,b,c,d);
	obj.BaselineStart=baseline(1,:);
	obj.BaselineLow=baseline(2,:);
	obj.BaselineUp=baseline(3,:);
end

% Update the goodness of fit
obj.RelStDev=gof.rmse;% due to the transformation to [0,1],
                      % the rmse is effectively a relative standard deviation
obj.CoeffDeterm=gof.rsquare;
obj.AdjCoeffDeterm=gof.adjrsquare;

% Update the performance stats of the fitting
obj.NumFunEvals=FitOutput.funcCount;
obj.NumIters=FitOutput.iterations;
obj.ExitFlag=FitOutput.exitflag;

end
