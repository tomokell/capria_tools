% This function displays (and outputs) an image in colormap gray,
% converting complex arrays into magnitude images, rotating if requested.
% The image is scaled to the percentiles of intensity MinPerc and MaxPerc.
% For 3D arrays, the image is displayed as a montage. For 4D arrays, the
% third dimension is concatenated along the horizontal direction and the
% fourth dimension along the vertical direction. NewWin creates a new
% figure window, FillWin makes the image as large as possible within the
% figure window, scales is a two component vector that allows images with
% anisotropic voxels to be displayed accurately. If UseAbsoluteScaling is
% true, MinPerc and MaxPerc are interpreted as absolute pixel values to
% scale the image before display, rather than as percentile values. If an
% output is requested, the scaled version of the image is returned, which
% can e.g. be saved to a file.
%
% Tom Okell, June 2022
%
% Im = DispIm(Im,MinPerc,MaxPerc,Rotate,NewWin,FillWin,scales,UseAbsoluteScaling)

function Im = DispIm(Im,MinPerc,MaxPerc,Rotate,NewWin,FillWin,scales,UseAbsoluteScaling)

if nargin < 2; MinPerc = 0; end
if nargin < 3; MaxPerc = 100; end
if nargin < 4; Rotate = false; end
if nargin < 5; NewWin = true; end
if nargin < 6; FillWin = false; end
if nargin < 7; scales = [1 1]; end
if nargin < 8; UseAbsoluteScaling = false; end

  % Open a new window
  if NewWin == true; figure; end
  
  % Rotate if necessary
  if Rotate == true; Im = rot90(Im); end
  
  % If complex, take the absolute value
  if ~isreal(Im); Im = abs(Im); end

  % If logical, convert to single
  if islogical(Im); Im = single(Im); end
  
  % Modify the image to allow display
  if UseAbsoluteScaling == true
      LowThr = MinPerc;
      UpThr = MaxPerc;
  else
      LowThr = prctile(Im(:),MinPerc);
      UpThr  = prctile(Im(:),MaxPerc);
  end
  
  if LowThr == UpThr; 
	  disp('DispIm: Warning - lower threshold = upper threshold')
	  disp('DispIm: Setting lower threshold to zero')
	  LowThr = 0; 
	  
	  if UpThr == 0;
		  disp('DispIm: Setting upper threshold to one')
		  UpThr = 1;
	  end
  end
  
  % Concatenate a 3D or 4D image (2D images are unaffected by this
  % function)
  if size(Im,4) == 1
      Im = CatSlices(Im,true); % Montage as well
  else
      Im = CatComps(Im);
  end
  
  % Create axis vectors
  x = (1:size(Im,2)) * scales(1);
  y = (1:size(Im,1)) * scales(2);
  
  % Display the image
  imagesc(x,y,Im,[LowThr, UpThr]);
  colormap gray;
  axis image; axis off;
  
  if FillWin == true
    % Make the image fill the window
    set(gca, 'Units', 'normalized', 'Position', [0 0 1 1]);
  end
  
  % Output the scaled version for saving
  if nargout > 0 
      Im(Im<LowThr) = LowThr;
      Im(Im>UpThr ) = UpThr;
      Im = (Im - LowThr)/(UpThr - LowThr);   
  end
  
  % Prevent spewing out numbers if no outputs are requested
  if nargout == 0; Im = []; end