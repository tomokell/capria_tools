% This function saves a complex array as magnitude and phase files using
% the Siemens phase convention (-pi to pi -> 0 to 4095)
%
% Tom Okell, June 2022
%
% function Save_Mag_Ph(ComplexIm, MagFName, PhFName, Scales)

function Save_Mag_Ph(ComplexIm, MagFName, PhFName, Scales)
  
  % Save the magnitude image
  disp('Saving magnitude image...')
  save_avw( abs(ComplexIm), MagFName, 'f', Scales );
  
  % Wrap the phase so that the interval -pi:0 becomes pi:2pi
  ph_im = angle(ComplexIm);
  ph_im( ph_im < 0 ) = ph_im( ph_im < 0 ) + 2*pi;
  
  % Save the phase image, converting from 0:2pi to 0:4095
  disp('Saving phase image...')
  save_avw( ph_im / (2 * pi) * 4095, PhFName, 'f', Scales );
  