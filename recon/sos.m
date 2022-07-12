% Perform a sqrt sum of squares of input x across dimension dim (assumed to be
% the final dimension of x if not specified).
%
% Mark Chiew, June 2022
%
% y = sos(x,dim)
function y = sos(x,dim)

N   =   ndims(x);
if nargin < 2
    dim = N;
end

y   =   sum(abs(x).^2,dim).^0.5;
y   =   permute(y, [setdiff(1:N, dim) dim]);
