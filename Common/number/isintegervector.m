function tf=isintegervector(x)
%% Is Input a Vector of Integer Numbers?
%
% See also: isintegerscalar, isintegermatrix, isintegerarray.
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 01/03/2013
% Last modified: 01/03/2013

tf=isintegerarray(x) && isvector(x);

end
