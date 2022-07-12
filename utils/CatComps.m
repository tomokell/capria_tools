% This function concatenates (four dimension) components and (third
% dimension) slices for display purposes
%
% Tom Okell, June 2022
%
% CatIm = CatComps(ImStack)
function CatIm = CatComps(ImStack)
  
  CatIm = [];
  for CompNo = 1:size(ImStack,4)
    CatIm = [CatIm; CatSlices(ImStack(:,:,:,CompNo))];
  end
