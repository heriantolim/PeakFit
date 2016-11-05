function p=transformpolycoeff(p,varargin)
%% Transform Polynomial Coefficient
%  This method transform a polynomial-coefficient vector p into another
%  polynomial-coefficient vector p' such that the corresponding polynomial
%  function y(x) transforms under: x->x'=ax+b and y->y'=cy+d.
%
% Syntax:
%  p=PeakFit.transformpolycoeff(p,a)
%  p=PeakFit.transformpolycoeff(p,a,b)
%  p=PeakFit.transformpolycoeff(p,a,b,c)
%  p=PeakFit.transformpolycoeff(p,a,b,c,d)
%
% Inputs
%  p: The polynomial coefficients in a row vector. It can be inputted as a 3-row
%     matrix, where the 1st, 2nd and 3rd row specifies the values of p, the
%     lower bound of p, and the upper bound of p respectively.
%
%  a, b, c, d: The linear transformation coefficiens.
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
% First created: 27/10/2016
% Last modified: 27/10/2016

assert(iscomplexmatrix(p) && (isrow(p) || size(p,1)==3),...
	'PeakFit:transformpolycoeff:InvalidInput',...
	['Input to the polynomial-coefficient vector must be a row vector or a ',...
		'3-row matrix of complex numbers.']);
assert(all(cellfun(@isrealscalar,varargin)),...
	'PeakFit:transformpolycoeff:InvalidInput',...
	'Input to the linear transformation coefficients must be a real scalar.');

N=nargin;
if N==1
	return
end

[m,n]=size(p);
a=varargin{1};
if a==0
	p=nan(m,n);
	return
end

if N>2
	b=varargin{2};
else
	b=0;
end

if N>3
	c=varargin{3};
	if c==0
		p=zeros(m,n);
		return
	end
else
	c=1;
end

if N>4
	d=varargin{4};
else
	d=0;
end

if n>1
	% compute the upper triangular matrix of the combinatorials
	A=zeros(n);
	N=n-1;
	for i=0:N
		for j=0:i
			A(n-i,n-i+j)=A(n-i,n-i+j)+nchoosek(i,j)*a^(-i)*(-b)^j;
		end
	end
	
	% transform the polynomial coefficients
	p(1,:)=p(1,:)*A;
	
	% transform the lower and upper bounds
	if m>2
		p2=zeros(1,n);
		p3=zeros(1,n);
		for j=1:n
			for i=1:j
				if A(i,j)<0
					p2(j)=p2(j)+p(3,j)*A(i,j);
					p3(j)=p3(j)+p(2,j)*A(i,j);
				else
					p2(j)=p2(j)+p(2,j)*A(i,j);
					p3(j)=p3(j)+p(3,j)*A(i,j);
				end
			end
		end
		p(2,:)=p2;
		p(3,:)=p3;
	end
end

p(:,n)=p(:,n)+d/c;
p=p*c;
if m>2 && c<0
	p([2,3],:)=p([3,2],:);
end

end