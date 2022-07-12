% Example script to reconstruct 4D CAPRIA angiographic images
%
% Tom Okell, June 2022

%% Check CAPRIA recon setup has been run first
if ~isvar('CAPRIA_recon_setup_complete')
    % Try running it if it's on the Matlab path
    disp('CAPRIA_recon_setup not yet run, trying to run now...')
    
    setupfname = which('CAPRIA_recon_setup');
    if isempty(setupfname)
        error('CAPRIA_recon_setup.m could not be found on the Matlab path')
    else
        run(setupfname)
    end
end

%% Set up parameters for the recon
MeasID = 284; % Assumes Siemens raw data in the format meas.dat format, with the measurement ID specified here
flags.regrid = false; % Don't do regridding
flags.calc_sens = true; % Calculate coil sensitivities as part of the reconstruction
flags.sense = true; % Perform a cgSENSE reconstruction
flags.L2 = false; % Don't perform an x sparse, L2 temporal smoothness recon
flags.xf = false; % Don't perform an xf sparse recon
flags.reconmean = false; % Don't reconstruct the mean image (we want the difference image here)
flags.LLR    = true; % Perform locally low rank recon
flags.xTV    = false; % Don't perform spatial total variation recon
flags.xTVtL2 = false; % Don't perform spatial total variation with L2 temporal smoothness recon
flags.LnS    = false; % Don't perform a low rank + sparse recon

SpatialMtxSz = 160; % Use the full spatial matrix size here to reconstruct at the maximum spatial resolution (1.1mm isotropic voxels)
TempResn = 216; % Temporal resolution for reconstruction (ms)
Nc = 12; % Number of coils to compress down to prior to reconstruction
L.S = 1E5*[1 1 1 1]; % Regularisation factor for cgSENSE
L.LLR.x =   1E-1; % Regularisation factor for LLR
L.LLR.p =   15; % LLR patch size

OutName = ['CFA_Angio']; % Output directory name

% Actual matrix size used for reconstruction, based on initial reconstructions to assess head sizes
Nx = round(SpatialMtxSz*230/160);
Ny = round(SpatialMtxSz*170/160);
Nz = round(SpatialMtxSz*248/160);
Nr = SpatialMtxSz*2; % Number of k-space points to include along the readout spoke (accounts for readout oversampling factor of 2)

% Define a shift of the centre of the image to account for the subject positioning
shift = floor([Nx Ny Nz]/2 + [10/230*Nx 0 -9/47*SpatialMtxSz]);

FrameNosForSensCalc = []; % Empty parameter means use all temporal frames for the coil sensitivity calculation
HannFilterFrac = [0.3]; % Apply heavier Hann filter which only spans the central 30% of k-space samples to minimize signal aliasing in this highly undersampled case.
SensThresh = [0.025]; % Threshold to mask coil sensitivities to exclude regions with no signal
SensKernelSize = []; % Use the default kernel size for adaptive combination coil sensitivity estimation
FracDataToUse = [1]; % Fraction of data to use for the reconstruction: set to one here to use all the data. Alternatively set to 0.5 to use the first half of the data, or [0.5 0.5] to use the second half of the data e.g. for repeatability analysis.

%% Run the reconstruction
CAPRIARecon('Raw_data',MeasID,'CAPRIA_Recon',OutName,false,TempResn,FracDataToUse,Nx,Ny,Nz,Nr,Nc,L,flags,SpatialMtxSz, ...
            shift,false,SensThresh,[],[],[],FrameNosForSensCalc,SensKernelSize,HannFilterFrac);

