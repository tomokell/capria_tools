function [psens, sens_xfm, s] = calc_psens(sens)
%   Mark Chiew
%   May 2014
%   Last Updated: Oct 2014
%
%   [psens, sens_xfm, s] = calc_psens(sens)
%
%   Calculate a compressed coil representation, psens, by deriving a linear
%   transform, sens_xfm, via an SVD of the original coil sensitivities,
%   sens. s is the normalised cumulative variance explained by the given
%   number of compressed coils.


dims            =   size(sens);
[u,s,sens_xfm]  =   lsvd(reshape(sens, [], dims(end)));
psens           =   reshape(u*s, dims);
s               =   cumsum(abs(diag(s)).^2)/sum(abs(diag(s)).^2);
