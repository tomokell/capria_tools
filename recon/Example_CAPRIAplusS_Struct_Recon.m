% Example script to reconstruct 4D CAPRIA+S structural images
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

%% Set up parameters for the recon - very similar to the angio recon
MeasID = 236; 
flags.regrid = false; 
flags.calc_sens = true; 
flags.sense = false; 
flags.L2 = false; 
flags.xf = false; 
flags.reconmean = true; % Reconstruct the mean signal to retain static tissue
flags.LLR    = true; % Only run LLR recon
flags.xTV    = false; 
flags.xTVtL2 = false; 
flags.LnS    = false;

SpatialMtxSz = 176;
TempResn = 327; % Use a larger temporal resolution than angio to better condition the recon
Nc = 8; 

% We use different regularization to the angio recon here
L.LLR.x =   1E-1;
L.LLR.p =   5;

OutName = [ns(MeasID) '_Struct']; % Output directory
CpName  = ['236_Angio']; % We can re-use the sensitivity maps from the angio recon

% Other parameters identical to the angio recon
Nx = round(SpatialMtxSz*84/64); 
Ny = SpatialMtxSz; 
Nz = round(SpatialMtxSz*56/64); 
Nr = SpatialMtxSz*2;
FrameNosForSensCalc = 5:6; % Not actually used here if we copy sensitivity maps from the angio recon
HannFilterFrac = [0.3]; 
SensThresh = [0.025];
SensKernelSize = [];

CAPRIARecon('.',MeasID,'../CAPRIA_Recon',OutName,false,TempResn,1,Nx,Ny,Nz,Nr,Nc,L,flags,SpatialMtxSz, ...
            [],false,SensThresh,[],[],[],FrameNosForSensCalc,SensKernelSize,HannFilterFrac,CpName);

