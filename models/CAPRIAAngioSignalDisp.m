% Returns the expected CAPRIA angiographic signal, accounting for
% dispersion, but ignoring the effect of RF pulses, M0b or inversion
% efficiency here.
%
% Tom Okell, May 2025
%
% Sig = CAPRIAAngioSignalDisp(t, v, delta_t, T1b, tau, s, p)
%
% where t is an array of time points to generate the expected signal, v is
% the fraction of the voxel occupied by macrovascular blood, delta_t is the
% macrovascular arrival time, T1b is the T1 of blood, tau is the labeling
% duration, s is the sharpness and p the time to peak of a gamma variate
% dispersion kernel. 

function Sig = CAPRIAAngioSignalDisp(t, v, delta_t, T1b, tau, s, p)

  % Calculate the modified parameters for the integral
  sprime = s + 1/T1b;
  
  % Calculate the scaling factor
  SF = 2 * v * exp(-delta_t/T1b) * (s/sprime)^(1+p*s);
  
  % Calculate the incomplete gamma integrals
  G = togammainc(sprime*(t-delta_t),1+p*s) - togammainc(sprime*(t-delta_t-tau),1+p*s);
  
  % Output the complete result
  Sig = SF * G;