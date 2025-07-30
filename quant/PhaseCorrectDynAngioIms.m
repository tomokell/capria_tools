% Calculates the mean phase of each voxel time series weighted by magnitude
% and corrects for it, outputting the real part only.
%
% RealCorrectedIms = PhaseCorrectDynAngioIms(Ims, TimeDim)
%
% Tom Okell, July 2025

function RealCorrectedIms = PhaseCorrectDynAngioIms(Ims, TimeDim)

  if nargin < 2; TimeDim = 3; end
  
  % Permute the images so time is the first dimension
  Dims = 1:6;
  PermOrder = [TimeDim Dims(Dims ~= TimeDim)];
  Ims = permute(Ims,PermOrder);
  
  % Initialise
  RealCorrectedIms = zeros(size(Ims));
  
  % Loop through the matrix, assuming time is the third dimension
  s = size(Ims);
  if length(s) < 6; s(6) = 1; end
  
  % Calculate the mean complex phase weighted by magnitude
  meanwphase = sum( Ims, 1 );
  meanwphase = meanwphase./abs(meanwphase);
  
  % Apply to the time series
  RealCorrectedIms = real(Ims .* repmat(conj(meanwphase),[s(1) 1 1 1 1 1]));
              
  % Perform the inverse permutation
  RealCorrectedIms = ipermute( RealCorrectedIms, PermOrder );