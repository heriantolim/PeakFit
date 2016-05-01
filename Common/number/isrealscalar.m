function tf=isrealscalar(x)
%% Is Input a Real Number (Scalar)?
% 
% See also: isrealvector, isrealmatrix, isrealarray.
% 
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 01/03/2013
% Last modified: 01/03/2013

tf=isrealarray(x) && isscalar(x);

end