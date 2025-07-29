% Returns the gamma inc function but first zeros in the elements of X which
% are negative
%
% Y = togammainc(X,A)
%
% Tom Okell, July 2025

function Y = togammainc(X,A,tail)

if nargin < 3; tail = 'lower'; end

X(X<0)=0; A(A<0) = 0;

Y = gammainc(X,A,tail);