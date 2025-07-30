% Overlay a mask on the specified image using the given colour
%
% OutIm = OverlayMask(Im,Mask,col)

function OutIm = OverlayMask(Im,Mask,col,NormaliseImFirst,ReplaceWithMask,NewWin)

if nargin < 3; col = [1 1 0]; end
if nargin < 4; NormaliseImFirst = true; end
if nargin < 5; ReplaceWithMask = false; end
if nargin < 6; NewWin = true; end

% Squeeze to ensure the third dimension is colour
Im = squeeze(Im);

% Normalise the image first if requested
if NormaliseImFirst
  Im = Im/max(Im(:));
end

% Check the third dimension has three colours
Sz = size(Im); if ndims(Im) < 3; Sz(3) = 1; end
if Sz(3) ~= 3
  if Sz(3) == 1  % Just pad to make a grayscale image
    Im = repmat(Im,[1 1 3]);
  else % Replace with VEPCASL colours
    Im = squeeze(Components2RGB(permute(Im,[1 2 4 3])));
  end
end

% Combine the image with the overlaid mask
ColMask = repmat(permute(col(:),[2 3 1]),Sz(1:2)).*repmat(Mask,[1 1 3]);

if ReplaceWithMask
  OutIm = Im;
  OutIm(repmat(Mask,[1 1 3])>0) = ColMask(repmat(Mask,[1 1 3])>0);  
else
  OutIm = Im + ColMask;
end

% Get rid of values above 1 and below 0
OutIm(OutIm>1) = 1;
OutIm(OutIm<0) = 0;

% Display
if NewWin; figure; end
imagesc(OutIm); axis equal; axis off;