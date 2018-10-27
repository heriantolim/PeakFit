function y=fnlorentzian(x,c,h,w)
%% Lorentzian function with center c, height h, and FWHM w
%
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 04/04/2013
% Last modified: 04/04/2013

y=h*(w/2)^2./((x-c).^2+(w/2)^2);

end