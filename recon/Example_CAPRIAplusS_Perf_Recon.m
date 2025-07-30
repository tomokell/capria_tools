% Example script to reconstruct 4D CAPRIA+S perfusion images
%
% Tom Okell, July 2025
%
% NB. This is very similar to the previous 4D CAPRIA recon script 
% "Example_CAPRIA_Perf_Recon.m", but with updated regularization 
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

SpatialMtxSz = 58; % Achieves an isotropic resolution of ~3.4mm
TempResn = 327; % Temporal resolution for reconstruction (ms)
Nc = 8; % Number of coils to compress down to prior to reconstruction
L.LLR.x =   1E-1; % Regularisation factor for LLR
L.LLR.p =   7; % LLR patch size

OutName = [ns(MeasID) '_Perf']; % Output directory name
CpName  = ['']; % If other reconstructions with the same FOV and matrix size have been run previously, you can avoid re-running e.g. 
                % the coil sensitivity calculation by specifying the directory name here

% Actual matrix size used for reconstruction, based on initial reconstructions to assess head sizes                
Nx = round(SpatialMtxSz*84/64); 
Ny = SpatialMtxSz; 
Nz = round(SpatialMtxSz*56/64); 
Nr = SpatialMtxSz*2; % Number of k-space points to include along the readout spoke (accounts for readout oversampling factor of 2)

FrameNosForSensCalc = 5:6; % Only use the last two frames for sensitivity calculation where the tissue magnetization is more stable
HannFilterFrac = []; % Use the default Hann filter
SensThresh = [0.025]; % Threshold to mask coil sensitivities to exclude regions with no signal
SensKernelSize = [3]; % Use a smaller kernel size for adaptive combination coil sensitivity estimation in this low resolution case

CAPRIARecon('.',MeasID,'../CAPRIA_Recon',OutName,false,TempResn,1,Nx,Ny,Nz,Nr,Nc,L,flags,SpatialMtxSz, ...
            [],false,SensThresh,[],[],[],FrameNosForSensCalc,SensKernelSize,HannFilterFrac,CpName);

