% Calculate the CAPRIA signal angiographic or perfusion signal with 
% a gamma variate dispersion kernel, assuming that all the blood 
% experiences all the RF pulses.
%
% Tom Okell, May 2025
%
% Usage: 
%   [dM, t, Alpha, R, dMTot, tAv, dMAv] = ...
%       CAPRIASignalDisp(AngioOrPerf,FAMode,FAParams,tstart,tau, ...
%             t0,Nsegs,Nphases,TR,physparams,CalcSigUsingMinAlpha, ...
%             OutputAsPerc,M0b,Alpha_inv)
%
%   Required inputs:
%       AngioOrPerf         =   'Angio' or 'Perf', as necessary 
%       FAMode              =   Flip angle mode (see CalcCAPRIAFAs)
%       FAParams            =   Flip angle parameters (see CalcCAPRIAFAs)
%       tstart              =   Start time of the simulation (s)
%       tau                 =   Labelling duration (s)
%       t0                  =   Start time of the readout (s)
%       Nsegs               =   Number of excitations within the temporal
%                               resolution (i.e. TempResn / TR)
%       Nphases             =   Number of frames (i.e. total readout time /
%                               (Nsegs*TR) )
%       TR                  =   Repetition time of the excitation pulses
%                               (s)
%       physparams          =   Physiological parameters struct. In all cases:
%                               .T1b = T1 of blood (s)
%                               For angiography:
%                               .v = macrovascular volume fraction
%                               .delta_t = macrovascular transit time (s)
%                               For perfusion:
%                               .f = CBF (ml/100g/min)
%                               .Deltat = arterial transit time to tissue (s)
%                               .T1 = T1 of tissue (s)
%                               For angio/perfusion with dispersion:
%                               .s = dispersion kernel sharpness (s^-1)
%                               .p = dispersion kernel time-to-peak (s)
%   Optional inputs:
%       CalcSigUsingMinAlpha=   Calculate the signal strength using the
%                               minimum flip angle within the readout
%                               (helps visualise the signal prior to the
%                               readout starting)
%       OutputAsPerc        =   Output signal strength as a percentage of
%                               M0b
%       M0b                 =   Equilibrium magnetisation of blood (au)
%       Alpha_inv           =   Inversion efficiency

function [dM, t, Alpha, R, dMTot, tAv, dMAv] = ...
            CAPRIASignalDisp(AngioOrPerf,FAMode,FAParams,tstart,tau, ...
            t0,Nsegs,Nphases,TR,physparams,CalcSigUsingMinAlpha, ...
            OutputAsPerc,M0b,Alpha_inv)

if nargin < 11; CalcSigUsingMinAlpha = false; end
if nargin < 12; OutputAsPerc = true; end
if nargin < 13; M0b = 1; end
if nargin < 14; Alpha_inv = 1; end

% Calculate a time array for the calculation. For simplicity, we separate
% the timepoints to be simulated by TR, which makes the RF attenuation
% calculation straightforward. First round the start time to make sure we
% can use multiples of TR and still hit t0 exactly.
tstart = t0 - round((t0 - tstart)/TR)*TR;

% Calculate the total readout time
Tro = TR*Nsegs*Nphases;

% Set up the time array
t = tstart:TR:(t0 + Tro - TR);

% Calculate flip angle scheme or scaling
Alpha = CalcCAPRIAFAs(FAMode,FAParams,t,t0);

% Work out whether to include dispersion or not
IncludeDisp = (isfield(physparams,'s') && isfield(physparams,'p'));

% Calculate the signal without attenuation
if strcmpi(AngioOrPerf,'angio') % Angiographic signal
    % Calculate the angio signal with no RF attenuation (no scaling for M0
    % or inversion efficiency)
    if ~IncludeDisp % No dispersion - simple version
        Sig = CAPRIAAngioSignalSimple(t, physparams.v, physparams.delta_t, physparams.T1b, tau);
    else % Include dispersion
        Sig = CAPRIAAngioSignalDisp  (t, physparams.v, physparams.delta_t, physparams.T1b, tau, physparams.s, physparams.p);
    end
    
elseif strcmpi(AngioOrPerf,'perf') % Perfusion signal
    % Calculate the Buxton model signal (no RF attenuation, no scaling for
    % M0 or inversion efficiency)
    if ~IncludeDisp % No dispersion - original Buxton model
        Sig =     BuxtonCASLModel(t,physparams.f,physparams.Deltat,tau,physparams.T1,physparams.T1b);
    else % Include dispersion
        Sig = BuxtonCASLModelDisp(t,physparams.f,physparams.Deltat,tau,physparams.s,physparams.p,physparams.T1,physparams.T1b);
    end
else % Unknown signal
    error('Unknown angio/perfusion type!')
end

% Scale the signal by M0b and the inversion efficiency
Sig = M0b * Alpha_inv * Sig;

% Calculate the RF attenuation factor
R = CAPRIAAttenuation(t,t0,Alpha);

% Adjust Alpha for signal calculation if necessary (creates smoother signal
% profiles for visualization)
AlphaSig = Alpha;
if CalcSigUsingMinAlpha
    AlphaSig(Alpha==0) = min(Alpha(Alpha>0));
end

% Calculate the final signal strength
dM = Sig .* R .* sin(AlphaSig*pi/180);

% Output as a percentage of M0 if requested
if OutputAsPerc
    dM = dM / M0b * 100; % Units of % of M0
end

% Calculate the summed measured signal during the readout if requested
if nargout > 4
    dMTot = sum(dM(t > t0));
end

% Average over output phases if requested
if nargout > 5
    tAv = zeros(Nphases,1);
    dMAv = zeros(Nphases,1);
    for ii = 1:Nphases
        Idx = (t >= t0+((ii-1)*Nsegs)*TR) & (t < t0+ii*Nsegs*TR);
        tAv(ii) = mean(t(Idx));
        dMAv(ii) = mean(dM(Idx));
    end
end
