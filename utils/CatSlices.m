% This function concatenates slices in a 3D volume to create an image for
% easier display.  If Montage is set to true, slices are concatenated in
% both directions rather than just horizontally. If Vertical is set to
% true, images are concatenated vertically rather than horizontally.
%
% Tom Okell, June 2022
%
% CatIm = CatSlices(InIm, Montage, Vertical)
function CatIm = CatSlices(InIm, Montage, Vertical)
  
  if nargin < 2; Montage = false; end
  if nargin < 3; Vertical = false; end
  
  % Determine the image size and pad with ones in case it's a low
  % dimensional image
  s = [size(InIm) 1 1];

  % Initialise output
  CatIm = [];

  % If not doing a montage, just concatenate slices in the
  % horizontal/vertical direction
  if Montage == false
    
    if Vertical
        for SliceNo = 1:size(InIm,3)
            CatIm = [CatIm; InIm(:,:,SliceNo,:,:)];
        end
    else        
        for SliceNo = 1:size(InIm,3)
            CatIm = [CatIm InIm(:,:,SliceNo,:,:)];
        end
    end
    
  else % Montage

    % Determine how many slices to concatenate in each direction
    NSlices = s(3);
    NImsX = ceil(sqrt(NSlices));
    NImsY = ceil(NSlices/NImsX);
    
    % Loop through the slices
    SlcNo = 1;
    for RowNo = 1:NImsY
      
      RowIm = [];
      
      for ColNo = 1:NImsX
        
        if SlcNo <= s(3)
          Slc = InIm(:,:,SlcNo,:,:);
        else
          Slc = zeros( [s(1:2) 1 s(4:5)] );
        end
        
        RowIm = [RowIm Slc];
        
        SlcNo = SlcNo + 1;
        
      end
      
      
      CatIm = [CatIm; RowIm];
        
    end
    
  end
