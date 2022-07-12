% Calculates a k-space trajectory for given azimuthal angles, Phis
% If Theta is also defined, a 3D trajectory will be produced. If NSamps
% is provided, it defines the number of k-space samples actually acquired
% for each line, which would be less than MtxSz for asymmetric echo.
% RawMtxSz is the full size of a full spoke, which may be different from
% the desired reconstructed matrix size. RawIdx refers to the indices which
% can be used to grab raw data corresponding to the sampled data points.
%
% Tom Okell, June 2022
%
% [kspace, RawIdx] = CalcRadialKSpaceTraj(Phi,MtxSz,Theta,NSamps,RawMtxSz,IsHalfSpoke)
function [kspace, RawIdx] = CalcRadialKSpaceTraj(Phi,MtxSz,Theta,NSamps,RawMtxSz,IsHalfSpoke)

if nargin < 3 || isempty(Theta); Theta = []; end
if nargin < 4 || isempty(NSamps);  NSamps = MtxSz; end
if nargin < 5 || isempty(RawMtxSz);  RawMtxSz = NSamps; end
if nargin < 6 || isempty(IsHalfSpoke);  IsHalfSpoke = false; end

% Set the number of spokes
Nspokes = length(Phi);

% Cope with partial Fourier/readout asymmetry/half-spoke radial here
if IsHalfSpoke
    kcentreIdx = 1;
    kCentreRecon = 1;
else
    kcentreIdx = ceil((RawMtxSz+1)/2); % Centre point of the full acquired spoke
    kCentreRecon = ceil((MtxSz+1)/2); % Centre point of the full spoke to be reconstructed
end
StartIdx = max(kcentreIdx - kCentreRecon, RawMtxSz - NSamps)+1;
EndIdx = kcentreIdx + (MtxSz - kCentreRecon);
ColIdx = (StartIdx:EndIdx)';
OutSamps = length(ColIdx);
if IsHalfSpoke
    krs = (ColIdx-kcentreIdx)/(MtxSz)*pi; % Radial k value for each spoke
else
    krs = (ColIdx-kcentreIdx)/(MtxSz/2)*pi; % Radial k value for each spoke
end

if isempty(Theta) % 2D

    % Initialise the trajectory
    kspace = zeros(Nspokes*OutSamps,2); % radians/pixel
    
    for ii = 1:Nspokes
        Idx = (ii-1)*OutSamps+1:ii*OutSamps;
        kspace(Idx,:) = krs * [sin(Phi(ii)) cos(Phi(ii))];
    end
    
else % 3D
    
    % Initialise the trajectory
    kspace = zeros(Nspokes*OutSamps,3); % radians/pixel
    
    for ii = 1:Nspokes
        Idx = (ii-1)*OutSamps+1:ii*OutSamps;
        kspace(Idx,:) = krs * [sin(Phi(ii))*sin(Theta(ii)) cos(Phi(ii))*sin(Theta(ii)) cos(Theta(ii))];
    end
    
end

% Output the index to be used on the raw k-space data
if nargout > 1
    % In the raw data the first sample is at point RawMtxSz - NSamps
    % so adjust indices accordingly
    RawIdx = ColIdx - (RawMtxSz - NSamps);
end