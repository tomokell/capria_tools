% Reads in a .nii.gz file (shorthand for read_avw).
%
% function [img, dims,scales,bpp,endian] = ra(fname)

function [img, dims,scales,bpp,endian] = ra(fname)
 
[img, dims,scales,bpp,endian] = read_avw(fname);
