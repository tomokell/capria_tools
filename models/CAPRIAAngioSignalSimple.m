% Very simple function to return the expected CAPRIA angiographic signal, not
% accounting for dispersion, the RF pulses, M0b or inversion efficiency here.
%
% Tom Okell, June 2022
%
% Sig = CAPRIAAngioSignalSimple(t, v, delta_t, T1b, tau)
%
% where t is an array of time points to generate the expected signal, v is
% the fraction of the voxel occupied by macrovascular blood, delta_t is the
% macrovascular arrival time, T1b is the T1 of blood and tau is the
% labeling duration.

function Sig = CAPRIAAngioSignalSimple(t, v, delta_t, T1b, tau)
    
    % Initialise output
    Sig = zeros(size(t));
    
    % Define the times at which the labeled blood is present in the voxel
    Idx = (t >= delta_t) & (t < delta_t + tau);
    
    % Calculate the signal
    Sig(Idx) = 2 * v * ones(size(Sig(Idx))) * exp(-delta_t / T1b);
end