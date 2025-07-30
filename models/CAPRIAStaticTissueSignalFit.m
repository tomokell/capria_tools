% Fit the CAPRIA+S static tissue signal to a T1 recovery model
%
% Tom Okell, July 2025
%
%   Usage:
%       [M0, InvAlpha, T1, B1rel, fit, Ims_pc] =
%           CAPRIAStaticTissueSignalFit(Ims,Mask,BGSMode,tau,BGStau1, ...
%           FAMode,FAParams,t0,TR,Nsegs,Nphases,IncB1rel,initmaps,...
%           fitmethod,CostRescale)
%
%   Outputs:
%       M0               =    Map of the fitted equilibrium magnetization
%       InvAlpha         =    Map of the fitted inversion pulse efficiency
%       T1               =    Map of the fitted T1
%       B1rel            =    Map of the fitted relative B1+ (if IncB1rel
%                             was set to true)
%       fit              =    Model fit to the data
%       Ims_pc           =    Phase corrected version of the input data
%                             used for fitting
%
%   Required Inputs:
%       Ims           =    4D complex CAPRIA+S tissue signal array
%       Mask          =    Only voxels within this 3D mask will be fit
%
%   Optional Inputs - parameters passed as empty will take the default values:
%       BGSMode       =    'Presat' for pre-saturation only or 'SI' for
%                           presat + a single inversion pulse (used in
%                           CAPRIA+S)
%       tau           =     PCASL label duration (s)
%       BGStau1       =     Time of the single inversion pulse (if used)
%                           relative to the start of PCASL labelling (s).
%       FAMode        =     Flip angle mode: 'CFA' (constant flip angle), 
%                           'Quadratic' or 'Maintain'
%       FAParams      =     Flip angle parameters (see CalcCAPRIAFAs.m)
%       t0            =     Time at which imaging starts relative to the 
%                           start of PCASL labelling (s)
%       TR            =     Time between excitation pulses in the readout
%                           (s)
%       Nsegs         =     Number of excitations within the temporal
%                           resolution (i.e. TempResn / TR)
%       Nphases       =     Number of frames (i.e. total readout time /
%                           (Nsegs*TR) )
%       IncB1rel      =     Include a relative B1+ scaling factor in the
%                           fitting process to accommodate B1+ inhomogeneity
%       initmaps      =     A structure containing maps to initialise the
%                           fitting: .M0, .InvAlpha, .T1 and (optionally) 
%                           .B1rel(ii,jj,kk)
%       fitmethod     =     The fitting algorithm passed to fmincon as
%                           opts.Algorithm
%       CostRescale   =     Rescaling approach for the cost function to
%                           ensure it does not get too small: 'global'
%                           rescales based on the maximum signal in the
%                           data, 'voxel' rescales each voxel separately or
%                           'none'. NB. Does not affect the scaling of the
%                           final output parameter maps.

function [M0, InvAlpha, T1, B1rel, fit, Ims_pc] = CAPRIAStaticTissueSignalFit(Ims,Mask,BGSMode,tau,BGStau1,FAMode,FAParams,t0,TR,Nsegs,Nphases,IncB1rel,initmaps,fitmethod,CostRescale)

if nargin < 3  || isempty(BGSMode); BGSMode = 'SI'; end
if nargin < 4  || isempty(tau); tau = 1.8; end
if nargin < 5  || isempty(BGStau1); BGStau1 = tau+2e-3+10e-3+11e-3/2; end
if nargin < 6  || isempty(FAMode); FAMode = 'quadratic'; end
if nargin < 7  || isempty(FAParams); FAParams = [2 9]; end
if nargin < 8  || isempty(t0); t0 = tau+2e-3+10e-3+11e-3+10e-3+2e-3; end
if nargin < 9  || isempty(TR); TR = 9.1e-3; end
if nargin < 10 || isempty(Nsegs); Nsegs = 36; end
if nargin < 11 || isempty(Nphases); Nphases = 6; end
if nargin < 12 || isempty(IncB1rel); IncB1rel = true; end
if nargin < 13 || isempty(initmaps); initmaps = []; end
if nargin < 14 || isempty(fitmethod); fitmethod = []; end
if nargin < 15 || isempty(CostRescale); CostRescale = 'voxel'; end

% Ensure the mask is logical and count voxels to fit
Mask = logical(Mask); Nvox = sum(Mask(:));

% Phase correct based on the final frame
Ims_pc = real(Ims.*exp(-1i*angle(Ims(:,:,:,end))));

% Define bounds for the fitting
if IncB1rel
    % params: [M0 InvAlpha T1  B1rel]
    lb =      [0    0.5   100e-3 0.9];
    ub =      [1e10 1.0   10.0   1.1];
else
    % params: [M0 InvAlpha  T1 ]
    lb =      [0    0.5  100e-3];
    ub =      [1e10 1.0    10.0];
end

% Initialise
fit = zeros(size(Ims_pc));
M0 = zeros(size(Ims_pc(:,:,:,1)));
InvAlpha = M0;
T1 = M0;
B1rel = M0;

% Constrained optimisation options
opts.Display = 'off';
if ~isempty(fitmethod)
    disp(['Using fit method: ' fitmethod])
    opts.Algorithm = fitmethod; % e.g. 'sqp'
end

% Loop through the data
VoxCount = 0;
for ii = 1:size(Ims_pc,1)
    for jj = 1:size(Ims_pc,2)
        for kk = 1:size(Ims_pc,3)
            
            if Mask(ii,jj,kk) == true
                disp([ns(VoxCount/Nvox*100) '% complete...'])
                   
                % Rescale the cost function if needed
                if strcmpi(CostRescale,'global')
                    % Global scale factor for fitting to ensure values don't get too small
                    scale = max(Ims_pc(:));
                elseif strcmpi(CostRescale,'voxel')
                    % Reset the scale for each voxel
                    scale = max(Ims_pc(ii,jj,kk,:));
                elseif strcmpi(CostRescale,'none')
                    % Don't rescale
                    scale = 1;
                else
                    error(['Unknown rescaling option: ' CostRescale]);
                end
                
                % Define the cost function: sum of squared differences
                % between data and model fit
                if IncB1rel
                    obfun = @(params) sum((squeeze(Ims_pc(ii,jj,kk,:))/scale - params(1)*CAPRIAStaticSignal(BGSMode,[BGStau1 params(2)],FAMode,FAParams*params(4),t0,params(3),TR,Nsegs,Nphases)/scale).^2 );                    
                else
                    obfun = @(params) sum((squeeze(Ims_pc(ii,jj,kk,:))/scale - params(1)*CAPRIAStaticSignal(BGSMode,[BGStau1 params(2)],FAMode,FAParams,t0,params(3),TR,Nsegs,Nphases)/scale).^2 );                    
                end
                
                % Initialise parameters
                if isempty(initmaps) % No maps provided: initialise with sensible parameters
                    if IncB1rel
                        params0 = [max(Ims_pc(ii,jj,kk,:)) 0.9 4.3 1.0];
                    else
                        params0 = [max(Ims_pc(ii,jj,kk,:)) 0.9 4.3];
                    end
                else % Use provided maps to intialise
                    if IncB1rel
                        params0 = [initmaps.M0(ii,jj,kk) initmaps.InvAlpha(ii,jj,kk) initmaps.T1(ii,jj,kk) initmaps.B1rel(ii,jj,kk)];
                    else
                        params0 = [initmaps.M0(ii,jj,kk) initmaps.InvAlpha(ii,jj,kk) initmaps.T1(ii,jj,kk)];
                    end
                    
                end

                % Run the constrained optimisation
                % Syntax: x = fmincon(fun,       x0, A, b,Aeq,beq,lb,ub,nonlcon,options)
                est         = fmincon(obfun,params0,[],[], [], [],lb,ub,     [],opts);
                
                % Store the outputs
                params = est;
                if ~IncB1rel; params(4) = 1; end

                % Model fit
                fit(ii,jj,kk,:) = params(1)*CAPRIAStaticSignal(BGSMode,[BGStau1 params(2)],FAMode,FAParams*params(4),t0,params(3),TR,Nsegs,Nphases);

                % Other parameters
                M0(ii,jj,kk) = params(1);
                InvAlpha(ii,jj,kk) = params(2);
                T1(ii,jj,kk) = params(3);
                B1rel(ii,jj,kk) = params(4);
                
                VoxCount = VoxCount + 1;
            end
        end
    end
end