function M = ifftdim(M,dim)

%
% [M] = ifftdim(m,dim)
%
%   Mark Chiew
%   May 2014
%
% performs centric Inverse Fourier Transform along dimension dim.
%

for i = dim
    M   =   fftshift(ifft(ifftshift(M, i), [], i), i)*sqrt(size(M,i));
end
