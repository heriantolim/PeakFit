function tf=isrealarray(x)
%% Is Input an Array of Real Numbers?
% 
% See also: isrealscalar, isrealvector, isrealmatrix.
% 
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 01/03/2013
% Last modified: 01/03/2013

tf=isnumeric(x) && isreal(x);

end