% Calculate flip angle schedules for CAPRIA acquisitions
% 
% Tom Okell, June 2022
%
% Usage:
%   Alpha = CalcCAPRIAFAs(FAMode,FAParams,t,t0)
%
% Required inputs:
%   FAMode      = 'CFA', 'Quadratic' or 'Maintain'
%   FAParams    = For CFA:          a scalar that defines the constant
%                                   flip angle in degrees.
%                 For Quadratic:    the flip angle varies quadratically
%                                   between FAParams(1) and FAParams(2) in
%                                   degrees.
%                 For Maintain:     Uses a backwards recursive formula to
%                                   maintain the signal at a constant level
%                                   (i.e. magnetisation loss in the
%                                   previous TR is counteracted by a higher
%                                   flip angle in the next TR). In this
%                                   case FAParams(1) defines the final flip
%                                   angle at the end of the readout.
%   t           = the time array to be simulated in s (assumes separation by TR)
%   t0          = the time at which imaging commences (s)

function Alpha = CalcCAPRIAFAs(FAMode,FAParams,t,t0)

% Initialise
Alpha = zeros(size(t));
Idx = t >=  t0;
N = sum(Idx); % Number of pulses played out

% CFA (FAParams = FA)
if strcmpi(FAMode,'CFA')
    Alpha(Idx) = FAParams(1); 
    
% VFA quadratic (FAParams = [FAMin FAMax])
elseif strcmpi(FAMode,'Quadratic')    
    Alpha(Idx) = FAParams(1) + (FAParams(2)-FAParams(1))*((0:(N-1))/(N-1)).^2;   
    
% VFA Maintain (FAParams = FAMax)
elseif strcmpi(FAMode,'Maintain')
    Alpha(Idx) = MaintainVFA(N,FAParams(1));
   
% Unknown
else
    error('Unknown FAMode!')
end