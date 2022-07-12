% Buxton (P)CASL model, from Buxton MRM 1998
%
% Tom Okell, June 2022
%
% dM = BuxtonCASLModel(t,f,deltat,tau,T1,T1b,lambda)
%
% Returns the difference in longitudinal magnetization for PCASL at
% timepoints after the start of labelling, t, given the CBF value, f
% (ml/100g/min), arterial transit time to tissue, Deltat, labelling
% duration, tau, T1 of tissue, T1, T1 of blood, T1b, and partition
% coefficient, lambda. tau can be given as a single value or one value for
% each element of t (for variable labelling duration experiments). Note
% that no account of M0b or the inversion efficiency is made here.

function dM = BuxtonCASLModel(t,f,Deltat,tau,T1,T1b,lambda)

if nargin < 5; T1 = 1.3; end
if nargin < 6; T1b = 1.6; end
if nargin < 7; lambda = 0.9; end

% Allow for single or variable tau
if length(tau) < length(t); tau = ones(size(t))*tau; end

% Scale the CBF parameter into s^-1
f = f/100/60;

% Calculate the effective T1 of tissue
T1p = 1/(1/T1+f/lambda);

% Calculate the approach to steady state term
qss = zeros(size(t));
Idx = (t>Deltat)&(t<=tau+Deltat);   qss(Idx) = 1-exp(-(t(Idx)-Deltat)/T1p);
Idx = t>tau+Deltat;                 qss(Idx) = 1-exp(-tau(Idx)/T1p);

% Calculate the longitudinal magnetisation difference
dM = zeros(size(t));
Idx = (t>=Deltat)&(t<Deltat+tau);   dM(Idx) = 2*f*T1p*exp(-Deltat/T1b)*qss(Idx);
Idx = (t>=Deltat+tau);              dM(Idx) = 2*f*T1p*exp(-Deltat/T1b)*exp(-(t(Idx)-tau(Idx)-Deltat)/T1p).*qss(Idx);

% Remove infinities and NaN values if they occur
dM(isinf(dM) | isnan(dM)) = 0;