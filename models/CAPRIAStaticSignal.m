% Calculate the CAPRIA static tissue signal (relative to M0)
%
% Tom Okell, July 2025
%
%   Usage: [Mxyav,Mxy,Mzminus] = CAPRIAStaticSignal(BGSMode,BGSParams,FAMode,FAParams,t0,T1,TR,Nsegs,Nphases)
%
%
%   Outputs:
%       Mxyav         =     Average transverse magnetization within each frame
%                           (phase)
%       Mxy           =     Transverse magnetization at each TR
%       Mzminus       =     Longitudinal magnetization just prior to each RF
%                           pulse
%
%   Required Inputs:
%       BGSMode       =    'Presat' for pre-saturation only or 'SI' for
%                           presat + a single inversion pulse (used in
%                           CAPRIA+S)
%       BGSParams     =     Only need for 'SI'. A two-element vector 
%                           defining the time of the single inversion pulse
%                           relative to the start of PCASL labelling (s)
%                           and the efficiency of the inversion pulse
%       FAMode        =     Flip angle mode: 'CFA' (constant flip angle), 
%                           'Quadratic' or 'Maintain'
%       FAParams      =     Flip angle parameters (see CalcCAPRIAFAs.m)
%       t0            =     Time at which imaging starts relative to the 
%                           start of PCASL labelling (s)
%       T1            =     Tissue T1 (s)
%       TR            =     Time between excitation pulses in the readout
%                           (s)
%       Nsegs         =     Number of excitations within the temporal
%                           resolution (i.e. TempResn / TR)
%       Nphases       =     Number of frames (i.e. total readout time /
%                           (Nsegs*TR) )

function [Mxyav,Mxy,Mzminus] = CAPRIAStaticSignal(BGSMode,BGSParams,FAMode,FAParams,t0,T1,TR,Nsegs,Nphases)

% Find the number of TRs
NTRs = Nsegs*Nphases;

% Initialise
Mzminus = zeros(NTRs,1); Mzplus = Mzminus; Mxy = Mzminus;

% Calculate flip angle scheme
Alpha = CalcCAPRIAFAs(FAMode,FAParams,(1:NTRs)*TR,TR);

% Start with zero Mz (assume pre-sat)
Mz0 = 0;

% Find the Mz at the start of the readout
switch lower(BGSMode)
    case 'presat' % Pure T1 recovery
        Mzminus(1) = T1recovery(Mz0,t0,T1);
        
    case 'si' 
        
        % Extra parameters from the input
        BGStau1 = BGSParams(1);
        BGSAlphaInv = BGSParams(2);
        
        % T1 recovery until the inversion pulse
        Mzpre180 = T1recovery(Mz0,BGStau1,T1);
        
        % Inversion
        Mzpost180 = Mzpre180*(1-2*BGSAlphaInv);
        
        % T1 recovery until the readout
        Mzminus(1) = T1recovery(Mzpost180,t0-BGStau1,T1);
        
    case 'di'
        error('DI not yet implemented');
    otherwise
        error('Unknown BGSMode')
end
        
% Loop through the imaging RF pulses
for ii = 1:NTRs
    if ii > 1 % Calculate the Mz- unless we're on the first pulse
        % T1 recovery since the last  pulse
        Mzminus(ii) = T1recovery(Mzplus(ii-1),TR,T1);
    end
    
    % RF attenuation of Mz
    Mzplus(ii) = Mzminus(ii)*cos(todeg2rad(Alpha(ii)));
    
    % Transverse magnetisation created (from Mzminus)
    Mxy(ii) = Mzminus(ii)*sin(todeg2rad(Alpha(ii)));
end

% Average across phases
Mxyav = zeros(Nphases,1);
for ii = 1:Nphases
    Idx = ((ii-1)*Nsegs)+1:ii*Nsegs;
    Mxyav(ii) = mean(Mxy(Idx));
end
