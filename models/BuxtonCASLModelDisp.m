% Buxton (P)CASL model, with dispersion, modified from Buxton MRM 1998
%
% Tom Okell, May 2025
%
% dM = BuxtonCASLModelDisp(t,f,deltat,tau,s,p,T1,T1b,lambda)
%
% Returns the difference in longitudinal magnetization for PCASL at
% timepoints after the start of labelling, t, given the CBF value, f
% (ml/100g/min), arterial transit time to tissue, deltat, labelling
% duration, tau, dispersion kernel sharpness, s (s^-1), dispersion kernel
% time-to-peak, p (s), T1 of tissue, T1 (s), T1 of blood, T1b (s), and
% partition coefficient, lambda. tau can be given as a single value or one
% value for each element of t (for variable labelling duration
% experiments). Note that no account of M0b or the inversion efficiency is
% made here.

function dM = BuxtonCASLModelDisp(t,f,Deltat,tau,s,p,T1,T1b,lambda)

if nargin < 7  || isempty(T1);  T1 = 1.3; end
if nargin < 8  || isempty(T1b); T1b = 1.6; end
if nargin < 9  || isempty(lambda); lambda = 0.9; end

% Assume M0b = 1 and inversion efficiency = 1
% Calculate some useful parameters
f = f/100/60; % Rescale CBF into s^-1 units
T1prime = 1/(1/T1+f/lambda); % Apparent tissue T1
sprime = s + 1/T1b; % Modified dispersion kernel sharpness
k = 1+p*s; % Gamma kernel parameter
beta = (1 - 1/(sprime*T1prime));

% Calculate the scaling factor
SF = 2 * f * (s/sprime)^k * T1prime * exp(-Deltat/T1b);
  
% Calculate the incomplete gamma integrals
Ga = togammainc(sprime*(t-Deltat),k) - togammainc(sprime*(t-Deltat-tau),k);
Gb = exp(-(t-Deltat)/T1prime)/beta^k.*(togammainc(beta*sprime*(t-Deltat),k) - exp(tau/T1prime)*togammainc(beta*sprime*(t-Deltat-tau),k));

% Final longitudinal magnetization difference
dM = SF .* (Ga - Gb);