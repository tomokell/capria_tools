% This function sets the nifti header of file Fname with the correct
% scanner orientation information using voxel sizes (mm) and temporal
% resolution (s), scales, image size, dims, and imaging region offset from
% isocentre, CentreSliceOffset (mm). Note that oblique imaging slices are
% not supported.
%
% Tom Okell, June 2022
%
% Mtx = tosetniftiorientinfo(Fname, scales, dims, CentreSliceOffset)
%
% Mtx converts voxel coords into x,y,z coords.

function Mtx = tosetniftiorientinfo(Fname, scales, dims, CentreSliceOffset)

  % Get the affine matrix
  Mtx = toaffineorientmatrix(scales,dims,CentreSliceOffset);
  
  % Convert to correct format
  Mtx = Mtx';
  Mtx = Mtx(:)';
  
  % Set the qform and sform matrices
  cmd = ['fslorient -setsformcode 1 ' Fname];
  disp(cmd); system(cmd);
  cmd = ['fslorient -setsform ' ns(Mtx) ' ' Fname];
  disp(cmd); system(cmd);
  cmd = ['fslorient -copysform2qform ' Fname];
  disp(cmd); system(cmd);