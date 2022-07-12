% Returns azimuthal and polar angles for N 3D golden means radial spokes as
% per Chan et al. MRM 2009.  If N is an array, this is taken as the golden
% ratio counter array. If ReverseOddLines is true, lines with an odd golden
% ratio counter are reversed (i.e. start at negative kz).
%
% Tom Okell, June 2022
%
% [Azi, Polar] = GoldenRatio3DAngles(N,ReverseOddLines)
%
function [Azi, Polar] = GoldenRatio3DAngles(N,ReverseOddLines)

if nargin < 2; ReverseOddLines = false; end

% Define the increments, from Chan et al
Phis = GoldenRatios3D;

% Calculate Polar and Azimuthal angles
% NB. apparent error in Chan paper figure here - Beta is the angle from the
% kz axis = Polar angle in Siemens terms
if length(N) == 1
    m = (0:(N-1))';
else
    m = N(:);
end

kz = mod(m*Phis(1),1); % Use GR to evenly distribute between 0 and 1

Polar = acos(kz); % Polar angle
Azi   = mod(m*Phis(2),1) * 2 * pi; % Azimuthal angle

% Reverse every other line if requested
if ReverseOddLines
    OddIdx = logical(mod(m,2));
    
    % Add pi to the azimuthal angle
    Azi(OddIdx) = mod( Azi(OddIdx) + pi, 2*pi);
    
    % Reverse kz
    Polar(OddIdx) = acos(-kz(OddIdx)); 
end
