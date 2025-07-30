% Example script to reconstruct 4D CAPRIA+S angiographic images
%
% Tom Okell, July 2025
%
% NB. This is very similar to the previous 4D CAPRIA recon script 
% "Example_CAPRIA_Angio_Recon.m", but with updated regularization 
% parameters and FOV in the slice direction due to the use of 
% slab-selective excitation.

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
MeasID = 236; % Assumes Siemens raw data in the format meas.dat format, with the measurement ID specified here
flags.regrid = false; % Don't do regridding
flags.calc_sens = true; % Calculate coil sensitivities as part of the reconstruction
flags.sense = false; % Don't perform a cgSENSE reconstruction
flags.L2 = false; % Don't perform an x sparse, L2 temporal smoothness recon
flags.xf = false; % Don't perform an xf sparse recon
flags.reconmean = false; % Don't reconstruct the mean image (we want the difference image here)
flags.LLR    = true; % Perform locally low rank recon
flags.xTV    = false; % Don't perform spatial total variation recon
flags.xTVtL2 = false; % Don't perform spatial total variation with L2 temporal smoothness recon
flags.LnS    = false; % Don't perform a low rank + sparse recon

SpatialMtxSz = 176; % Use the full spatial matrix size here to reconstruct at the maximum spatial resolution (1.1mm isotropic voxels)
TempResn = 163; % Temporal resolution for reconstruction (ms)
Nc = 8;  % Number of coils to compress down to prior to reconstruction
L.LLR.x =   7E-2; % Regularisation factor for LLR
L.LLR.p =   15; % LLR patch size

OutName = [ns(MeasID) '_Angio']; % Output directory name
CpName  = ['']; % If other reconstructions with the same FOV and matrix size have been run previously, you can avoid re-running e.g. 
                % the coil sensitivity calculation by specifying the directory name here

% Actual matrix size used for reconstruction, based on initial reconstructions to assess head sizes
Nx = round(SpatialMtxSz*84/64); 
Ny = SpatialMtxSz; 
Nz = round(SpatialMtxSz*56/64); % Note this is a lot smaller than previous 4D CAPRIA recons due to the use of a slab-select excitation pulse
Nr = SpatialMtxSz*2; % Number of k-space points to include along the readout spoke (accounts for readout oversampling factor of 2)

FrameNosForSensCalc = 10:12; % Only use the last few frames for sensitivity calculation where the tissue magnetization is more stable
HannFilterFrac = [0.3]; % Apply heavier Hann filter which only spans the central 30% of k-space samples to minimize signal aliasing in this highly undersampled case.
SensThresh = [0.025]; % Threshold to mask coil sensitivities to exclude regions with no signal
SensKernelSize = []; % Use the default kernel size for adaptive combination coil sensitivity estimation

CAPRIARecon('.',MeasID,'../CAPRIA_Recon',OutName,false,TempResn,1,Nx,Ny,Nz,Nr,Nc,L,flags,SpatialMtxSz, ...
            [],false,SensThresh,[],[],[],FrameNosForSensCalc,SensKernelSize,HannFilterFrac,CpName);

