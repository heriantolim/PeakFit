function y=fngaussian(x,c,h,w)
%% Gaussian function with center c, height h, and FWHM w
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 04/04/2013
% Last modified: 04/04/2013

y=h*exp(-4*log(2)*(x-c).^2/w^2);

end