function [B,index]=sortcolumns(A,varargin)
%% Sort the Columns of a Matrix
%  B=sortcolumns(A) sorts the columns of A in ascending order of the first row.
%  Strings are sorted in the familiar dicitionary order, which uses ASCII
%  values. When the first column has equal values, sortcolumns sorts according
%  to the next row and repeats this behavior for succeeding equal values.
%
%  [B,index]=sortcolumns(A) additionally returns a vector of the permutation
%  of indices.
%
%  This function takes extra arguments in the similar way as the sortrows
%  function.
%
% See also: sortrows.
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 24/03/2013
% Last modified: 25/03/2013

[B,index]=sortrows(A',varargin{:});
B=B';

end