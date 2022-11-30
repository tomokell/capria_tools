%   Combined angiography and perfusion using radial imaging and ASL (CAPRIA)
%   reconstruction code, built on Mark Chiew's operators and the IRT NUFFT
%
%   Tom Okell (tokell@fmrib.ox.ac.uk) and Mark Chiew (mchiew@fmrib.ox.ac.uk)
%   June 2022
%
%   NB: Requires the nufft portion of Jeff Fessler's irt toolbox
%       See http://web.eecs.umich.edu/~fessler/irt/fessler.tgz
%       Saving out data in Nifti format also requires some FSL utilities: 
%       See https://fsl.fmrib.ox.ac.uk/fsl/
%       Siemens raw data reading relies on mapVBVD, by Philipp Ehses and
%       colleagues.
%
%   Usage:
%       CAPRIARecon(RawDir,MeasID,OutDir,OutName,Seq2D,TempResn,FracDataToUse, ...
%           Nx,Ny,Nz,Nr,Nc,L,flags,SpatialMtxSz,shift,DispResults,SensitivityThresh, ...
%           SeqType,PhaseCorr,FrameNos,FrameNosForSensCalc,SensKernelSize, ...
%           HannFilterFrac,SensCpDir,PreCompressCoils,IsHalfSpoke)   
%
%   Required Inputs:
%       RawDir        =    Directory which contains the raw meas.dat file
%       MeasID        =    Numeric measurement ID of the meas.dat file
%       OutDir        =    Directory to save output images
%       OutName       =    Subdirectory name to create to save output images
%       Seq2D         =    Is the sequence 2D? Leave empty to determine
%                          this from the raw data file
%       TempResn      =    Temporal resolution (ms) to reconstruct at
%       FracDataToUse =    Fraction of data to use in the reconstruction.
%                          Specify as a single value (e.g. 0.5 would use
%                          the first half of the raw data) or a two
%                          component vector [Frac Offset], to specify the
%                          fraction of data to use and an offset from the
%                          beginning of the file (e.g. [0.5 0.5] would use
%                          half the data, starting from halfway through).
%       Nx            =    Number of voxels along x in the reconstructed image
%       Ny            =    Number of voxels along y in the reconstructed image
%       Nz            =    Number of voxels along z in the reconstructed image
%       Nr            =    Number of readout points along each spoke to
%                          include in the reconstruction (including readout 
%                          oversampling)
%       Nc            =    Number of virtual coils to compress down to
%                          prior to running the reconstruction
%
%   Optional Inputs - parameters passed as empty will take the default values:
%       L                   =   Regularization factors for reconstruction,
%                               specified as a struct, e.g. for cgSENSE: 
%                               L.S = 1E4*[1 1 1 1]; % Components are x, y, z and t
%                               For LLR: L.LLR.x = 1E-2, L.LLR.p = 5, etc.
%                               Any missing parameters take the default
%                               values.
%       flags               =   Binary flags to determine which (parts of)
%                               reconstructions to perform, specified as a
%                               struct:
%                               flags.regrid: perform regridding (adjoint
%                                             NUFFT)
%                               flags.calc_anat: Calculate "anatomical"
%                                                image, using mean label/
%                                                control data. Used for
%                                                coil sensitivity
%                                                estimation
%                               flags.calc_sens: Calculate coil
%                                                sensitivities
%                               flags.sense:     perform cgSENSE recon
%                               flags.L2:        perform x-sparse, L2 temporal
%                                                smoothness recon
%                               flags.xf:        perform xf-sparse recon
%                               flags.reconmean: recon only the mean image
%                                                (i.e. average label and control 
%                                                data prior to recon rather 
%                                                than taking the difference)
%                               flags.recontag:  recon only the tag data
%                               flags.reconcontrol: recon only the control
%                                                   data
%                               flags.LLR:       perform locally low rank
%                                                recon
%                               flags.xTV:       perform spatial total
%                                                variation recon
%                               flags.xTVtL2:    perform spatial TV with L2
%                                                temporal smoothness recon
%                               flags.LnS:       perform low rank + sparse
%                                                recon
%       SpatialMtxSz        =   nominal matrix size to reconstruct at
%                               (determines reconstructed voxels size via
%                               nominal FOV (from raw data) / SpatialMtxSz
%       shift               =   FOV shift (default is 50% image size)
%       DispResults         =   Display results during recon (set to false
%                               for cluster usage)
%       SensitivityThresh   =   Threshold to mask coil sensitivity maps
%       SeqType             =   CAPRIA sequence type: default is to read
%                               from the raw data header, which works for
%                               all current CAPRIA sequences.
%       PhaseCorr           =   perform phase correction step to align
%                               k-space data acquired in label/control
%                               conditions
%       FrameNos            =   Frame numbers to reconstruct
%       FrameNosForSensCalc =   Frame numbers to use for coil sensitivity
%                               calculations
%       SensKernelSize      =   Kernel size for adaptive combine coil
%                               sensitivity estimation
%       HannFilterFrac      =   Fraction of the k-space spoke to apply the
%                               Hann filter to (zeros elsewhere) for
%                               anatomical image calculation (used in coil
%                               sensitivity estimation).
%       SensCpDir           =   If a previous reconstruction has already
%                               estimated coil sensitivities, these can be
%                               copied from the specified directory to prevent 
%                               recalculating them
%       PreCompressCoils    =   Option to compress coils based on an SVD of the
%                               raw k-space data, prior to coil sensitivity
%                               estimation
%       IsHalfSpoke         =   For radial-out (half spoke) data (limited
%                               testing - use with care)
%       
%
function CAPRIARecon(RawDir,MeasID,OutDir,OutName,Seq2D,TempResn,FracDataToUse,Nx,Ny,Nz,Nr,Nc,L,flags,SpatialMtxSz,shift,DispResults,SensitivityThresh,SeqType,PhaseCorr,FrameNos,FrameNosForSensCalc,SensKernelSize,HannFilterFrac,SensCpDir,PreCompressCoils,IsHalfSpoke)

% Deal with optional parameters
if nargin < 13 || isempty(L)
    L = [];
end
if nargin < 14 || isempty(flags)
    flags = [];
end
if nargin < 15
    SpatialMtxSz = [];
end
if nargin < 16 || isempty(shift)
    shift = [];
end
if nargin < 17 || isempty(DispResults)
    DispResults = false;
end
if nargin < 18 || isempty(SensitivityThresh)
    SensitivityThresh = 0.1;
end
if nargin < 19 || isempty(SeqType)
    SeqType = [];
end
if nargin < 20 || isempty(PhaseCorr)
    PhaseCorr = true;
end
if nargin < 21 || isempty(FrameNos)
    FrameNos = [];
end
if nargin < 22 || isempty(FrameNosForSensCalc)
    FrameNosForSensCalc = [];
end
if nargin < 23 || isempty(SensKernelSize)
    SensKernelSize = 5;
end
if nargin < 24 || isempty(HannFilterFrac)
    HannFilterFrac = 1;
end
if nargin < 25 || isempty(SensCpDir)
    SensCpDir = [];
end
if nargin < 26 || isempty(PreCompressCoils)
    PreCompressCoils = false;
end
if nargin < 27 || isempty(IsHalfSpoke)
    IsHalfSpoke = false;
end

%% Deal with optional parameters in fields

% Regularisation parameters
% Specify default values here
Ldef.S = 1E5*[1 1 1 1]; 
Ldef.L2.x  =   1E-2;
Ldef.L2.t  =   1E4;
Ldef.xf    =   1E-4;
Ldef.LLR.x =   1E-1;
Ldef.LLR.p =   15;
Ldef.xTV   =   1E-3;
Ldef.xTVtL2.x = 1E-3;
Ldef.xTVtL2.t = 1E4;
Ldef.LnS.L = 5E-9;
Ldef.LnS.S = 2E-10;

% Assign user-specified parameters, if available, otherwise resort to the
% default values
L = topassargs(L,Ldef,true,true);

% Repeat for flags
flagsdef.regrid = true;
flagsdef.calc_anat = true;
flagsdef.calc_sens = true;
flagsdef.sense  = true;
flagsdef.L2     = false;
flagsdef.xf     = false;
flagsdef.reconmean = false;
flagsdef.recontag = false;
flagsdef.reconcontrol = false;
flagsdef.LLR    = true;
flagsdef.xTV    = false;
flagsdef.xTVtL2 = false;
flagsdef.LnS    = false;

flags = topassargs(flags,flagsdef,true,true);


    
%% Create output directories
if ~exist(OutDir)
    disp(['Output directory ' OutDir ' does not exist: creating...'])
    mkdir(OutDir);
end
OutSubDir = [OutDir '/' OutName];
if ~exist(OutSubDir)
    disp(['Output directory ' OutSubDir ' does not exist: creating...'])
    mkdir(OutSubDir);
end

%% Display the inputs for the log file and save for reference
ReconStartDateTime = datetime; 
disp(ReconStartDateTime)
disp('>>> Starting CAPRIARecon')
disp('> Using the following parameters:')
disp(['RawDir = ' ns(RawDir)])
disp(['MeasID = ' ns(MeasID)])
disp(['OutDir = ' ns(OutDir)])
disp(['OutName = ' ns(OutName)])
disp(['Seq2D = ' ns(Seq2D)])
disp(['TempResn = ' ns(TempResn)])
disp(['FracDataToUse = ' ns(FracDataToUse)])
disp(['Nx = ' ns(Nx)])
disp(['Ny = ' ns(Ny)])
disp(['Nz = ' ns(Nz)])
disp(['Nr = ' ns(Nr)])
disp(['Nc = ' ns(Nc)])
disp(['L = '])
disp(L)
disp(['flags = '])
disp(flags)
disp(['SpatialMtxSz = ' ns(SpatialMtxSz)])
disp(['shift = ' ns(shift)])
disp(['DispResults = ' ns(DispResults)])
disp(['SensitivityThresh = ' ns(SensitivityThresh)])
disp(['SeqType = ' ns(SeqType)])
disp(['PhaseCorr = ' ns(PhaseCorr)])
disp(['FrameNos = ' ns(FrameNos)])
disp(['FrameNosForSensCalc = ' ns(FrameNosForSensCalc)])
disp(['SensKernelSize = ' ns(SensKernelSize)])
disp(['HannFilterFrac = ' ns(HannFilterFrac)])
disp(['SensCpDir = ' SensCpDir])
disp(['PreCompressCoils = ' ns(PreCompressCoils)])
disp(['IsHalfSpoke = ' ns(IsHalfSpoke)])
save([OutSubDir '/params.mat']);

%% Copy anat and sens data from another directory if requested
CurDir = pwd;
if ~isempty(SensCpDir)
    disp(['Creating symbolic links to previous anat and sens data inside ' SensCpDir '...'])
    cd(OutSubDir)
    
    % Check the copy directory exists
    if ~exist(['../' SensCpDir])
        error(['Could not find directory to copy data from: ' SensCpDir])
    end
    
    % Double check that the relevant parameters used for this recon match
    % those used for the previous sensitivity calculation (image size
    % and shift)
    if exist(['../' SensCpDir '/params.mat'])
       tmp = load('params.mat');
       newparams = [Nx,Ny,Nz,Nr,SpatialMtxSz,shift];
       cpparams = [tmp.Nx,tmp.Ny,tmp.Nz,tmp.Nr,tmp.SpatialMtxSz,tmp.shift];
       if sum(newparams ~= cpparams) > 0
           disp(['Current critical parameters: ' newparams])
           disp(['Copy critical parameters: ' cpparams])
           error('Critical parameters do not match - cannot copy sens/anat data - aborting')
       end
    else
       warning('Sens copy directory does not contain params.mat - assuming parameters match...') 
    end
    
    % Create symbolic links to avoid using more disk space than is
    % necessary
    tosystem(['ln -s ../' SensCpDir '/anat.mat'])
    tosystem(['ln -s ../' SensCpDir '/sens.mat'])
    
    % Go back to original directory
    cd(CurDir)
end

%% Read in data
disp('Reading in data...')

% Get into the raw data directory
cd(RawDir)

% Extract the relevant data for this reconstruction (including separate
% label and control data)
tic
[kdata,k,S,CentreSlcOffset,scales] = ExtractCAPRIAData(MeasID,Seq2D,TempResn,FracDataToUse,false,SpatialMtxSz,SeqType,PhaseCorr,FrameNos,[],IsHalfSpoke);
% Dimensions of kdata are now: kspace_samples time Coils label/control
% Dimensions of k are now:     kspace_samples time xyz
% S contains the dimensions of the raw data
% CentreslcOffset is the physical displacement of the centre of the FOV
% from isocentre
% scales are the voxel dimensions and temporal resolution

% Take the mean of label and control data
kdata_mean = mean(kdata,4);

% Deal with options to reconstruct only the mean data, control data or
% label data, otherwise take the difference between label and control
if flags.reconmean
    kdata_diff = kdata_mean;
elseif flags.recontag
    kdata_diff = kdata(:,:,:,1);
elseif flags.reconcontrol
    kdata_diff = kdata(:,:,:,2);
else
    kdata_diff = kdata(:,:,:,2)-kdata(:,:,:,1);
end

% Clear data to save memory
clear kdata;

% Report timing
t = toc; disp(['-- Data reading took ' ns(t/60/60) ' hours'])


%% Define useful parameters for later based on raw data
NcOrig = size(kdata_mean,3);
Nt = size(kdata_mean,2);
if Seq2D
    Ndims = 2;
else
    Ndims = 3;
end

% If not specified, use all frames for coil sensitivity calculations
if isempty(FrameNosForSensCalc)
    FrameNosForSensCalc = 1:Nt;
end

%% Pre-compress coils if requested based on the raw k-space data directly
if PreCompressCoils
    tic
    disp('Pre-compressing coils...')
    
    % Calculate the transformation to compress the coils
    [~, xfm, ~] = calc_psens(reshape(kdata_mean(:,FrameNosForSensCalc,:),[],NcOrig));
    
    % Apply the coil compression to the difference and mean data
    kdata_diff = apply_sens_xfm(xfm, kdata_diff,Nc,3);
    kdata_mean = apply_sens_xfm(xfm, kdata_mean,Nc,3);
    
    % Reset the number of coils for later processing
    NcOrig = Nc;
    
    % Report the timing
    t = toc; disp(['-- Coil pre-compression took ' ns(t/60/60) ' hours'])
end

%% Calculate a smoothed anatomical image from mean data
cd(CurDir); cd(OutSubDir);

% If it already exists, then load in
if exist('./anat.mat')
    tic
    disp('Loading anat data...')
    load ./anat.mat
    t = toc; disp(['-- Anat loading took ' ns(t/60/60) ' hours'])
    
% Otherwise, generate if it's been requested    
elseif flags.calc_anat
    tic
    disp('Constructing anatomical image...')
    
    % Define the NUFFT operator
    if isempty(shift)
        E = xfm_NUFFT([Nx, Ny, Nz, 1],[],[],reshape(k(:,FrameNosForSensCalc,:),[],1,Ndims),'wi',1,'table',true); 
    else
        E = xfm_NUFFT([Nx, Ny, Nz, 1],[],[],reshape(k(:,FrameNosForSensCalc,:),[],1,Ndims),'wi',1,'shift',shift,'table',true);
    end
    
    % Apply a Hann filter to smooth the data - no high res info needed for
    % coil sensitivity estimation
    Nf = round(Nr*HannFilterFrac); % Shrink the filter if requested
    if mod(Nf,2) == 0; Nf = Nf - 1; end % Make sure it's odd, so we can centre about k = 0
    kcentreIdx = ceil((Nr+1)/2); % Find the k=0 index
    fIdx = (kcentreIdx - floor(Nf/2)):(kcentreIdx + floor(Nf/2)); % Calculate the index centred about this point
    
    % Calculate the filter
    F = zeros(Nr,1);
    F(fIdx) = hann(Nf);
    
    % Modify to cope with partial Fourier
    if S(1) < Nr % Asymmetric echo was used, so just keep the end part of the filter
        F = F((Nr-S(1)+1):end);
    end
    
    % Apply the filter to the mean data
    kdata_mean = reshape(reshape(kdata_mean(:,FrameNosForSensCalc,:),length(F),[],NcOrig).*F,[],NcOrig);
    
    % Regrid each coil separately using an iterative procedure
    for i = 1:NcOrig
        disp(['Regridding coil ' ns(i)])
        m(:,i) = E.iter(kdata_mean(:,i),@pcg,1E-4,50,[1E3,1E3,1E3,0]);
    end
    
    % Reshape the output to the image size x number of coils
    m = reshape(m,[E.Nd NcOrig]);
    
    % Calculate a sum of squares image
    SoSM = sos(m);
    
    % Display if requested
    if DispResults
        DispIm(squeeze(flip(SoSM,1))); title 'Anat'
    end
    
    % Save the output to a .mat file
    q = matfile('anat','Writable',true);
    q.m = m;
    q.anat = SoSM;
    
    % Save to Nifti
    SaveCAPRIAToNifti(SoSM,'anat',scales,CentreSlcOffset,[],shift);
    
    % Report timing
    t = toc; disp(['-- Anatomical construction took ' ns(t/60/60) ' hours'])
    
else % Do not perform anatomical image reconstruction (not needed for simple regridding, for example)
    disp('Skipping anatomical image construction...')
end

% Save memory
clear E kdata_mean

%% Estimate Sensitivities and Compress
% If they already exist, load them in
if exist('./sens.mat')
    tic
    disp('Loading sensitivity data...')
    load ./sens.mat
    t = toc; disp(['-- Sens loading took ' ns(t/60/60) ' hours'])
    
% Otherwise, calculate if requested
elseif flags.calc_sens
    tic
    disp('Estimating coil sensitivities...')
    
    % Run the adaptive combine algorithm
    sens = adaptive_estimate_sens('data',m,'kernel',SensKernelSize,'thresh',SensitivityThresh,'verbose',true);
    
    % Save the output to a .mat file
    save('sens','sens','-v7.3');
    
    % ...and a Nifti file
    SaveCAPRIAToNifti(sens,'sens',scales,CentreSlcOffset,[],shift);
    
    % Report timing
    t = toc; disp(['-- Coil sensitivity estimation took ' ns(t/60/60) ' hours'])
    
else % Don't calculate sensitivities
    disp('Skipping coil sensitivity calculation...')
    sens = [];
end

% Calculate and apply coil compression
if isempty(sens) % No sensitivity calculation, so we just leave the data as it is here
    ps = []; % Leave compressed coil variable empty
    pd = kdata_diff; % Just copy the kspace data without compression
    Nc = size(pd,3); % Reset the number of coils to be the full set of physical coils
    
else % We have coil information, so compress (if needed) and apply
    
    if ~PreCompressCoils % If we haven't used pre-compression, compress here and apply
        tic
        disp('Calculating compressed coils...')
        [ps, xfm, ~] = calc_psens(sens); % Creates ps, the compressed coils
    
        disp('Applying coil sensitivites...')
        pd = apply_sens_xfm(xfm, kdata_diff,Nc,3); % Creates pd, the coil-compressed k-space data
        
        t = toc; disp(['-- Coil compression and application took ' ns(t/60/60) ' hours'])
        
    else % Otherwise, just copy the data, which is already compressed
        disp('No coil compression needed at this stage as performed previously');
        ps = sens;
        pd = kdata_diff;
    end
        
end

% Save memory
clear m sens kdata_diff xfm

%% Regridding for comparison
% If already performed, then don't repeat
if exist('./resr.mat')
    disp('Regridding recon already done...')

% Otherwise, regrid if requested
elseif flags.regrid
    tic
    disp('Performing regridding reconstruction...')
    
    % Define the NUFFT operator
    if isempty(shift)
        E0 = xfm_NUFFT([Nx,Ny,Nz,Nt],[], [], k,'table',true);
    else
        E0 = xfm_NUFFT([Nx,Ny,Nz,Nt],[], [], k,'shift',shift,'table',true);
    end
    
    % Regrid each coil using the adjoint operator. Note that the
    % weights calculated by the operator (E0.w) are defined as
    % the square root of the density compensation weights, and are applied
    % in both forward and adjoint transforms. Therefore, for measured
    % k-space data we need to pre-multiply by E0.w, so that when the
    % adjoint operator is applied, it includes an addition factor of E0.w,
    % thereby giving the full density compensation required.
    for i = 1:Nc
        disp(['Regridding coil ' ns(i) ' of ' ns(Nc)])
        imgr(:,:,:,i) = E0'*(pd(:,:,i).*E0.w);
    end
    
    % Reshape to the correct dimensions
    imgr = reshape(imgr,[E0.Nd E0.Nt Nc]);
    
    % Coil combination
    if isempty(ps) % No coil information: Use sum of squares
        imgr = sos(imgr);
    else % We have coil information, so do Roemer combination    
        imgr = sum(imgr.*conj(reshape(ps(:,:,:,1:Nc),[E0.Nd 1 Nc])),5)./sum(ps(:,:,:,1:Nc).*conj(ps(:,:,:,1:Nc)),4);
    end
    
    % Save the regridding result as a .mat file
    save('resr','imgr','-v7.3')
    
    % ...and to Nifti
    SaveCAPRIAToNifti(abs(imgr),'resr',scales,CentreSlcOffset,[],shift);
    
    % Report timing
    t = toc; disp(['-- Regridding recon took ' ns(t/60/60) ' hours'])
    
    % Display if requested
    if DispResults
        DispIm(mean(flip(imgr(:,:,1:ceil(end/10):end,:),1),4),1,99.5); title 'Regridding temporal mean'
        DispIm(flip(squeeze(imgr(:,:,1:ceil(end/10):end,1:ceil(end/4):end)),1),1,99.5); title 'Regridding'
    end

end

% Save memory
clear E0 imgr

%% Pre-calculate the encoding operator, PSF etc. memory savings
% Only calculate if running a more advanced recon below
if flags.sense || flags.L2 || flags.xf || flags.LLR || flags.xTV || flags.xTVtL2 || flags.LnS
    tic
    disp('Pre-calculations for memory savings...')
    
    % If already calculated in a previous run, the reload here
    if exist('./data.mat')
        disp('Loading data...')
        load ./data.mat
        
        % Recreate the encoding operator, E, from saved data (this is faster than recreating from scratch):
        if isempty(shift)
            E = xfm_NUFFT([Nx,Ny,Nz,Nt],ps(:,:,:,1:Nc), [], k, 'wi',1,'table',true, 'st', [], 'PSF', PSF);
        else
            E = xfm_NUFFT([Nx,Ny,Nz,Nt],ps(:,:,:,1:Nc), [], k, 'wi',1,'shift',shift,'table',true, 'st', [], 'PSF', PSF);
        end
        
    % Otherwise, generate from scratch
    else
       
        % Generate the encoding operator
        if isempty(shift)
            E = xfm_NUFFT([Nx,Ny,Nz,Nt],ps(:,:,:,1:Nc), [], k, 'wi',1,'table',true);
        else
            E = xfm_NUFFT([Nx,Ny,Nz,Nt],ps(:,:,:,1:Nc), [], k, 'wi',1,'shift',shift,'table',true);
        end
        
        % Calculate the adjoint transform applied to the coil compressed k-space
        % data
        dd = E'*pd(:,:,1:Nc);

        % Save dd and the PSF to a .mat file to save time if re-running
        % recons later
        q = matfile('./data','Writable',true);
        
        q.dd = dd;
        q.PSF = E.PSF;
    end
    
    % Report timing
    t = toc; disp(['-- Precalculations took ' ns(t/60/60) ' hours'])
end

% Clear unnecessary variables - these are all encapsulated inside E now
clear ps pd k PSF


%% cgSENSE reconstruction with regularisation
% Only run if not already performed
if exist('./resS.mat')
    disp('SENSE recon already done...')

% Otherwise, run if requested
elseif flags.sense
    disp('Performing SENSE recon...')

    % Hardwire the number of iterations here
    Nits = 100;
    
    % Time how long one iteration takes, then use this to predict the full
    % computational time (usually overestimates by quite a bit)
    tic; 
    E.iter(dd, @pcg, 1E-4, 1, L.S); 
    OneIt = toc; 
    disp(['One iteration took ' ns(OneIt) ' s. Estimated time for ' ns(Nits) ' iterations is ' ns(OneIt*Nits) ' s = ' ns(OneIt*Nits/60/60) ' hrs'])
    
    % Run the main iterative recon
    tic
    imgS = E.iter(dd, @pcg, 1E-4, Nits, L.S);
    
    % Reshape the output
    imgS = reshape(imgS,[E.Nd E.Nt]);
    
    % Save to .mat
    save('resS','imgS','-v7.3')
    
    % Save to Nifti
    SaveCAPRIAToNifti(abs(imgS),'resS',scales,CentreSlcOffset,[],shift);

    % Report timing
    t = toc; disp(['-- SENSE recon took ' ns(t/60/60) ' hours'])
    
    % Display if requested
    if DispResults
        DispIm(mean(flip(imgS(:,:,1:ceil(end/10):end,:),1),4),1,99.5); title 'SENSE temporal mean'
        DispIm(flip(squeeze(imgS(:,:,1:ceil(end/10):end,1:ceil(end/4):end)),1),1,99.5); title 'SENSE'
    end

end

% Save memory
clear imgS


%% x-sparse + temporal L2 smoothness recon
% Check if the recon has already been run
if exist('./resL2.mat')
    disp('L2 recon already done...')
    
% Otherwise, run if requested
elseif flags.L2
    tic
    disp('Performing L2 recon...')

    % Hardwire the iterations here
    iter    =   200;
    
    % Calculate the maximum step size and divide by two to be on the safe
    % side
    disp('Calculating step size...')
    step = E.max_step()/2;

    % Run the recon
    imgL2 =   fista_L2(E, dd, L.L2.x, L.L2.t, [prod(E.Nd) E.Nt], iter, step);

    % Reshape the output
    imgL2 =   reshape(imgL2,[E.Nd E.Nt]);
    
    % Save as a .mat
    save('resL2','imgL2','-v7.3')
    
    % Save as Nifti
    SaveCAPRIAToNifti(abs(imgL2),'resL2',scales,CentreSlcOffset,[],shift);

    % Report timing
    t = toc; disp(['-- L2 recon took ' ns(t/60/60) ' hours'])

    % Display if requested
    if DispResults
        DispIm(mean(flip(imgL2(:,:,1:ceil(end/10):end,:),1),4),1,99.5); title 'L2 temporal mean'
        DispIm(flip(squeeze(imgL2(:,:,1:ceil(end/10):end,1:ceil(end/4):end)),1),1,99.5); title 'L2'        
    end

end

% Save memory
clear imgL2

%% xf-sparse recon
% Check if recon was run already
if exist('./resxf.mat')
    disp('xf recon already done...')

% Otherwise, run if requested
elseif flags.xf
    tic
    disp('Performing xf recon...')

    % Hardwire iterations
    iter    =   200;
    
    % Calculate step size, if it wasn't done already
    if ~exist('step','var')
        disp('Calculating step size...')
        step = E.max_step()/2;  
    end
    
    % Run the recon
    imgxf =   fista_xf(E, dd, L.xf, [prod(E.Nd) E.Nt], iter, step);
    
    % Reshape the output
    imgxf =   reshape(imgxf,[E.Nd E.Nt]);
    
    % Save to .mat
    save('resxf','imgxf','-v7.3');
    
    % Save to Nifti
    SaveCAPRIAToNifti(abs(imgxf),'resxf',scales,CentreSlcOffset,[],shift);
    
    % Report timing
    t = toc; disp(['-- xf recon took ' ns(t/60/60) ' hours'])

    % Display if requested
    if DispResults
        DispIm(mean(flip(imgxf(:,:,1:ceil(end/10):end,:),1),4),1,99.5); title 'xf temporal mean'
        DispIm(flip(squeeze(imgxf(:,:,1:ceil(end/10):end,1:ceil(end/4):end)),1),1,99.5); title 'xf'        
    end
    
end

%% Locally Low Rank (LLR) recon
% Check if already run
if exist('./resLLR.mat')
    disp('LLR recon already done...')

% Otherwise, run the recon if requested    
elseif flags.LLR
    tic
    disp('Performing LLR recon...')
    
    % Hardwire the iterations
    iter    =   200;
    
    % Run the recon
    imgLLR =   pogm_LLR(E, dd, L.LLR.x, [1 1 1]*L.LLR.p, [E.Nd E.Nt], iter);
    
    % Reshape the output
    imgLLR =   reshape(imgLLR,[E.Nd E.Nt]);
    
    % Save to .mat
    save('resLLR','imgLLR','-v7.3');
    
    % Save to Nifti
    SaveCAPRIAToNifti(abs(imgLLR),'resLLR',scales,CentreSlcOffset,[],shift);
    
    % Report timing
    t = toc; disp(['-- LLR recon took ' ns(t/60/60) ' hours'])

    % Display if requested
    if DispResults
        DispIm(mean(flip(imgLLR(:,:,1:ceil(end/10):end,:),1),4),1,99.5); title 'LLR temporal mean'
        DispIm(flip(squeeze(imgLLR(:,:,1:ceil(end/10):end,1:ceil(end/4):end)),1),1,99.5); title 'LLR'        
    end
    
end

%% Spatial total variation (xTV) recon
% Check if recon was already run
if exist('./resxTV.mat')
    disp('xTV recon already done...')
    
% Otherwise, run the recon if requested    
elseif flags.xTV
    tic
    disp('Performing xTV recon...')
    
    % Hardwire the iterations
    iter    =   100;
    
    % Run the recon
    imgxTV =   fgp_xTV(E, dd, L.xTV, iter);
    
    % Reshape the output
    imgxTV =   reshape(imgxTV,[E.Nd E.Nt]);
    
    % Save to .mat
    save('resxTV','imgxTV','-v7.3');
    
    % Save to Nifti
    SaveCAPRIAToNifti(abs(imgxTV),'resxTV',scales,CentreSlcOffset,[],shift);
    
    % Report timing
    t = toc; disp(['-- xTV recon took ' ns(t/60/60) ' hours'])

    % Display if requested
    if DispResults
        DispIm(mean(flip(imgxTV(:,:,1:ceil(end/10):end,:),1),4),1,99.5); title 'xTV temporal mean'
        DispIm(flip(squeeze(imgxTV(:,:,1:ceil(end/10):end,1:ceil(end/4):end)),1),1,99.5); title 'xTV'        
    end

end


%% Spatial total variation (xTV) + L2 over time recon
% Check if recon was already run
if exist('./resxTVtL2.mat')
    disp('xTVtL2 recon already done...')
    
% Otherwise, run the recon if requested    
elseif flags.xTVtL2
    tic
    disp('Performing xTVtL2 recon...')
    
    % Hardwire the iterations
    iter    =   100;
    
    % Run the recon
    imgxTVtL2 =   fgp_xTV_tL2(E, dd, L.xTVtL2.x, L.xTVtL2.t, iter);
    
    % Reshape the output
    imgxTVtL2 =   reshape(imgxTVtL2,[E.Nd E.Nt]);
    
    % Save to .mat
    save('resxTVtL2','imgxTVtL2','-v7.3');
    
    % Save to Nifti
    SaveCAPRIAToNifti(abs(imgxTVtL2),'resxTVtL2',scales,CentreSlcOffset,[],shift);

    % Report timing
    t = toc; disp(['-- xTVtL2 recon took ' ns(t/60/60) ' hours'])

    % Display if requested
    if DispResults
        DispIm(mean(flip(imgxTVtL2(:,:,1:ceil(end/10):end,:),1),4),1,99.5); title 'xTVtL2 temporal mean'
        DispIm(flip(squeeze(imgxTVtL2(:,:,1:ceil(end/10):end,1:ceil(end/4):end)),1),1,99.5); title 'xTVtL2'        
    end
    
end

%% Low rank and sparse recon
% Check if recon was already run
if exist('./resLnS.mat')
    disp('L+S recon already done...')
    
% Otherwise, run the recon if requested    
elseif flags.LnS
    tic
    disp('Performing L+S recon...')
    
    % Hardwire the iterations
    iter    =   100;
    
    % Run the recon
    [imgLnS_L,imgLnS_S] =   pogm_LS(E, dd, L.LnS.L, L.LnS.S, [E.Nd E.Nt], iter);
        
    % Combine the low rank (L) and sparse (S) parts and reshape the output
    imgLnS =   reshape(imgLnS_L + imgLnS_S,[E.Nd E.Nt]);
    
    % Clear the separate components to save memory
    clear imgLnS_L imgLnS_S
    
    % Save to .mat
    save('resLnS','imgLnS','-v7.3');
    
    % Save to Nifti
    SaveCAPRIAToNifti(abs(imgLnS),'resLnS',scales,CentreSlcOffset,[],shift);
    
    % Report timing
    t = toc; disp(['-- L+S recon took ' ns(t/60/60) ' hours'])

    % Display if requested
    if DispResults
        DispIm(mean(flip(imgLnS(:,:,1:ceil(end/10):end,:),1),4),1,99.5); title 'LnS temporal mean'
        DispIm(flip(squeeze(imgLnS(:,:,1:ceil(end/10):end,1:ceil(end/4):end)),1),1,99.5); title 'LnS'        
    end
    
end



%% Finish!
% Go back to original directory
cd(CurDir)

disp('DONE!')
