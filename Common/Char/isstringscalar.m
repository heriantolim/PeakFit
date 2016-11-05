function tf=isstringscalar(x)
%% Is Input a String?
% 
% See also: isstringvector, isstringmatrix, isstringarray.
% 
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 01/03/2013
% Last modified: 04/03/2013

tf=ischar(x) && (isrow(x) || isempty(x));

end