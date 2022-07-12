% Calculate the attenuation due to RF pulses in a CAPRIA-style acquisition.
%
% Tom Okell, June 2022
%
% R = CAPRIAAttenuation(t,t0,Alpha)
%
% where t is an array of timepoints separated by the TR, t0 is the start of
% imaging and Alpha is an array of flip angles of size(t).

function R = CAPRIAAttenuation(t,t0,Alpha)

% Initialise
R = zeros(size(t)); R(1) = 1; 

% Calculate attenuation due to each previous RF pulse
for ii = 2:length(t)
    if t(ii) > t0
        R(ii) = R(ii-1)*cos(todeg2rad(Alpha(ii-1))); % Attenuation
    else
        R(ii) = 1;
    end
end