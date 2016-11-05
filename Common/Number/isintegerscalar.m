function tf=isintegerscalar(x)
%% Is Input an Integer Number (Scalar)?
% 
% See also: isintegervector, isintegermatrix, isintegerarray.
% 
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 01/03/2013
% Last modified: 01/03/2013

tf=isintegerarray(x) && isscalar(x);

end