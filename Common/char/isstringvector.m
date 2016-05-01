function tf=isstringvector(x)
%% Is Input a Cell Array of Strings with Size 1-by-N or M-by-1?
% 
% See also: isstringscalar, isstringmatrix, isstringarray.
% 
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 01/03/2013
% Last modified: 01/03/2013

tf=isstringarray(x) && isvector(x);

end