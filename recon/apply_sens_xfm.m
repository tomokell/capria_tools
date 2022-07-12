function dd = apply_sens_xfm(xfm, dd, p, coil_dim)
%   
%   transformed_data = apply_sens_xfm(xfm, data, p, coil_dim)
%
%   Mark Chiew
%   May 2014
%
%   Apply a coil compression transform, xfm, to the input data, compressing
%   down to p coils. coil_dim is the coil dimension of the input data.

dims    =   size(dd);
dims(3) =   size(dd,3);
dims(4) =   size(dd,4);
dims(coil_dim) = p;

x       =   [setdiff(1:4, coil_dim) coil_dim];
y(x)    =   1:4;
dd      =   permute(reshape(reshape(permute(dd,x),[],size(dd,coil_dim))*xfm(:,1:p),dims(x)),y);
