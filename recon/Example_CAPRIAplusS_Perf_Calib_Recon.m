% Example script to reconstruct 4D CAPRIA+S perfusion calibration images
% by reconstructing the mean signal rather than the control-label difference
% signal, but with a matrix size matched to the perfusion recon.
%
% Tom Okell, July 2025

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
MeasID = 236; 
flags.regrid = false; 
flags.calc_sens = true; 
flags.sense = false; 
flags.L2 = false; 
flags.xf = false; 
flags.reconmean = true; % Reconstruct the mean signal, not the difference
flags.LLR    = true; % Only use LLR recon
flags.xTV    = false; 
flags.xTVtL2 = false; 
flags.LnS    = false;

% Match the perfusion recon parameters
SpatialMtxSz = 58;
TempResn = 327;
Nc = 8; 
L.LLR.x =   1E-1;
L.LLR.p =   7;

OutName = [ns(MeasID) '_Perf_Calib'];

% Copy the sensitivity maps from the previously run standard perf recon to avoid
% unnecessary identical processing
CpName  = ['236_Perf'];

% Match perfusion recon parameters
Nx = round(SpatialMtxSz*84/64); 
Ny = SpatialMtxSz; 
Nz = round(SpatialMtxSz*56/64); 
Nr = SpatialMtxSz*2;

FrameNosForSensCalc = 5:6; 
HannFilterFrac = []; 
SensThresh = [0.025];
SensKernelSize = [3];

CAPRIARecon('.',MeasID,'../CAPRIA_Recon',OutName,false,TempResn,1,Nx,Ny,Nz,Nr,Nc,L,flags,SpatialMtxSz, ...
            [],false,SensThresh,[],[],[],FrameNosForSensCalc,SensKernelSize,HannFilterFrac,CpName);

