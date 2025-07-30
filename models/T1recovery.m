% T1 recovery for a time period delta, starting at Mz0
%
% Tom Okell, July 2025
%
% Mz = T1recovery(Mz0,delta,T1)
%
% Where Mz0 is the starting longitudinal magnetization, delta is the time
% step and T1 is the longitudinal relaxation time

function Mz = T1recovery(Mz0,delta,T1)

Mz = 1 - (1 - Mz0).*exp(-delta/T1);