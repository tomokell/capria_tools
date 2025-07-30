% Calibrate reconstructed CAPRIA data to allow perfusion to be calculated
% in absolute units. Since CAPRIA isn't fully integrated into the oxasl
% pipeline yet, we assume here there is a reference output directory (e.g.
% from comparison 3D-GRASE data) which contains useful masks etc. produced
% by oxasl with CSF calibration, but this could also be reproduced by
% comparison with a T1.anat directory from a separate T1w structural or
% CAPRIA+S structural image.
%
% Tom Okell, July 2025
%
%   Usage: CAPRIACalibration(PerfDir,PerfCalibDir,RefDir,OutDir, ...
%                            UseRefBias,ReconFName,CalibReconFName,...
%                            BGSMode,tau,BGStau1,FAMode,FAParams,t0,...
%                            TR,Nsegs,Nphases,TE)
%
%   Required Inputs:
%       PerfDir         =   CAPRIA reconstruction directory containing the
%                           perfusion images.
%       PerfCalibDir    =   CAPRIA reconstruction directory containing the
%                           calibration images to be used (identical to the
%                           perfusion recon but the mean signal rather than
%                           the difference signal reconstructed)
%       RefDir          =   Reference oxasl output directory containing
%                           useful masks etc.
%       OutDir          =   Output directory
%
%   Optional inputs:
%       UseRefBias      =   Use the bias field from the reference directory
%                           to correct the CAPRIA data
%       ReconFName      =   Name of the reconstructed images in the CAPRIA
%                           recon directory (e.g. 'resLLR')
%       CalibReconFName =   Name of the reconstructed images in the CAPRIA
%                           calibration recon directory (e.g. 'resLLR')
%       BGSMode         =   'Presat' for pre-saturation only or 'SI' for
%                           presat + a single inversion pulse (used in
%                           CAPRIA+S)
%       tau             =   Labelling duration (s)
%       BGStau1         =   Only need for 'SI'. The time of the single
%                           inversion pulse relative to the start of PCASL
%                           labelling (s)
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
%       TE            =     Echo time (s)

function CAPRIACalibration(PerfDir,PerfCalibDir,RefDir,OutDir,UseRefBias,ReconFName,CalibReconFName,BGSMode,tau,BGStau1,FAMode,FAParams,t0,TR,Nsegs,Nphases,TE)

if nargin < 5; UseRefBias = false; end
if nargin < 6; ReconFName = 'resS'; end
if nargin < 7; CalibReconFName = 'resS'; end
if nargin < 8; BGSMode = 'SI'; end
if nargin < 9; tau = 1.8;  end
if nargin < 10; BGStau1 = tau+2e-3+10e-3+11e-3/2;  end
if nargin < 11; FAMode = 'quadratic';  end
if nargin < 12; FAParams = [2 9];  end
if nargin < 13; t0 = tau+2e-3+10e-3+11e-3+10e-3; end
if nargin < 14; TR = 9.1e-3;  end
if nargin < 15; Nsegs = 36;  end
if nargin < 16; Nphases = 6; end
if nargin < 17; TE = 4.74e-3; end

% Remove any extensions to the reconstructed file name
ReconFName      = regexprep(regexprep(     ReconFName,'\.nii\.gz',''),'\.mat',''); 
CalibReconFName = regexprep(regexprep(CalibReconFName,'\.nii\.gz',''),'\.mat',''); 

% Make the output directory
mkdir(OutDir)

% Convert inputs to absolute paths to avoid complications with relative
% path names
PerfDir = char(getAbsPath(PerfDir));
PerfCalibDir = char(getAbsPath(PerfCalibDir));
RefDir = char(getAbsPath(RefDir));
OutDir = char(getAbsPath(OutDir));
CurDir = pwd;

% Move into the output directory
cd(OutDir)

% Create symbolic links to the input data
tosystem(['ln -s ' PerfCalibDir '/' CalibReconFName '.nii.gz calib_mag.nii.gz'])
tosystem(['ln -s ' PerfCalibDir '/' CalibReconFName '.mat calib.mat'])
tosystem(['ln -s ' PerfDir      '/' ReconFName '.mat perf.mat'])

%% Register the ventricle mask from oxasl to CAPRIA space
if ~exist('ventriclemask.nii.gz')
    tosystem(['fslmaths calib_mag.nii.gz -Tmean calib_mag_Tmean'])
    tosystem(['flirt -dof 6 -usesqform -in ' RefDir '/calibration/calib_img.nii.gz -ref calib_mag_Tmean.nii.gz -applyxfm -out ref2CAPRIA -omat ref2CAPRIA.mat'])
    tosystem(['flirt -in ' RefDir '/calibration/refmask.nii.gz -ref calib_mag_Tmean.nii.gz -init ref2CAPRIA.mat -applyxfm -out ventricles'])
    tosystem(['fslmaths ventricles.nii.gz -thr 0.95 -bin ventriclemask.nii.gz'])
end

%% Repeat for GM and WM - have to transform PVE estimates
if ~exist('gm_mask.nii.gz') || ~exist('wm_mask.nii.gz')
    tosystem(['convert_xfm -omat struc2CAPRIA.mat -concat ref2CAPRIA.mat ' RefDir '/reg/struc2asl.mat'])
    tosystem(['to_applywarp ' RefDir '/structural/gm_pv.nii.gz calib_mag_Tmean.nii.gz gm_pv_aslspace struc2CAPRIA.mat'])
    tosystem(['to_applywarp ' RefDir '/structural/wm_pv.nii.gz calib_mag_Tmean.nii.gz wm_pv_aslspace struc2CAPRIA.mat'])
    tosystem(['fslmaths gm_pv_aslspace.nii.gz -thr 0.5 -bin gm_mask'])
    tosystem(['fslmaths wm_pv_aslspace.nii.gz -thr 0.9 -bin wm_mask'])
end

%% Bias correction
if UseRefBias % Register the reference directory bias estimate into CAPRIA space
    tosystem(['flirt -in ' RefDir '/senscorr/sensitivity.nii.gz -ref calib_mag_Tmean.nii.gz -init ref2CAPRIA.mat -applyxfm -out bias']);
    
else % Run fsl_anat on the CAPRIA structural low res data to get a bias field
    % Need to remove NaN values first or FSL tools fail. Also rescale to be on
    % the safe side.
    tosystem(['fslmaths calib_mag_Tmean.nii.gz -nan -mul 1e9 calib_mag_Tmean_rescaled'])

    % Run fsl_anat
    tosystem(['fsl_anat --nocrop --nosubcortseg --noreorient -i calib_mag_Tmean_rescaled.nii.gz'])

    % NB. the bias field here appears to only be the refinement stage, so
    % doesn't correct for the main bias in the image, so calculate manually instead
    tosystem(['fslmaths calib_mag_Tmean_rescaled.anat/T1.nii.gz -div calib_mag_Tmean_rescaled.anat/T1_biascorr.nii.gz calib_mag_Tmean_rescaled_bias'])
end

%% Load the calibration images
calib = load('calib.mat'); calibfnames = fieldnames(calib);
eval(['calib = calib.' calibfnames{1} ';'])

%% Bias correct
if UseRefBias
    bias = ra('bias');
else
    bias = ra('calib_mag_Tmean_rescaled_bias');
end

% Fill in the zeros
bias(bias==0) = 1;

% Reorient to CAPRIA convention
bias = flipdim(bias,3);
bias = flipdim(bias,2);
PermOrder = [2 1 3 4 5];
bias = ipermute(bias,PermOrder);

% Divide through by the bias field
calib_bc = calib./repmat(bias,[1 1 1 size(calib,4)]);

% Show the pre and post- bias correction in some slices
slcs = round(linspace(10,size(calib,3),6));
DispIm([calib(:,:,slcs,end); calib_bc(:,:,slcs,end)])
title 'Original calibration images (above), bias corrected (below)'

%% Grab the masks and fit to the model for each tissue type
masknames = {'ventriclemask','gm_mask','wm_mask'};
LargeFigWindow(1,1);

for jj = 1:length(masknames)    
    mask = logical(ra(masknames{jj}));

    % Reorient to CAPRIA convention
    mask = flipdim(mask,3);
    mask = flipdim(mask,2);
    PermOrder = [2 1 3 4 5];
    mask = ipermute(mask,PermOrder);
    
    % Find slices which contain some mask
    MaskSlcIdx = find(squeeze(sum(sum(mask)))>0);
    slcs = round(linspace(min(MaskSlcIdx),max(MaskSlcIdx),6));

    % Display an overlay of the mask on the calibration image
    subplot(2,3,jj);    
    OverlayMask(CatSlices(abs(calib_bc(:,:,slcs,end)),true),CatSlices(mask(:,:,slcs),true),[1 1 0],true,false,false);

    % Phase correct 
    calib_bc_pc = real(calib_bc.*exp(-1i*angle(calib_bc(:,:,:,end))));    

    % Extract signals from within this mask at each timepoint
    MeanSigPC = [];
    for ii = 1:size(calib_bc_pc,4)
        tmp = calib_bc_pc(:,:,:,ii);
        MeanSigPC(ii) = mean(tmp(mask));
    end    
    
    % Define the timing
    TempResn = Nsegs*TR;
    t = ((1:size(calib_bc_pc,4))-0.5)*TempResn;

    % Define the cost function for fitting (sum of squared differences
    % between model and data)
    % NB. We scale up in the cost function to avoid early stopping of the
    % fit due to the small mean value of the data, but this doesn't affect
    % the final scaling parameter output (M0)
    obfun = @(params) sum((MeanSigPC(:)*1e10 - params(1)*CAPRIAStaticSignal(BGSMode,[BGStau1 params(2)],FAMode,FAParams,t0,params(3),TR,Nsegs,Nphases)*1e10).^2 );
    
    % Bounds and initialisation
    lb = [0    0.5 100e-3];
    ub = [1e10 1.0 10.0];
    params0 = [max(MeanSigPC) 0.9 4.3];
    
    % Run the constrained optimisation
    % Syntax: x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
    est = fmincon(obfun,params0,[],[],[],[],lb,ub,[],[]);
    
    % Calculate the model fit
    params = est;
    fit = params(1)*CAPRIAStaticSignal(BGSMode,[BGStau1 params(2)],FAMode,FAParams,t0,params(3),TR,Nsegs,Nphases);

    % Display
    subplot(2,3,jj+3)
    plot(t,MeanSigPC,'o',t,fit,'-');
    M0(jj) = est(1);
    InvAlpha(jj) = est(2);
    T1(jj) = est(3);
    xlabel 'Time/s'
    ylabel 'Mean signal/au'
    legend('Data','Fit')
    title([masknames{jj} ': M0 = ' num2str(M0(jj),3) ', T1 = ' num2str(T1(jj),3) ', InvAlpha = ' num2str(InvAlpha(jj),3)],'interpreter','none')
end

%% Calculate M0b using WM as a reference - larger signal, more voxels
jj = find(contains(masknames,'wm'));
if isempty(jj)
    jj = find(contains(masknames,'WM'));
end

% T2* correction, as per FSL's asl_calib:
T2starref = 50e-3; PC = 0.82; % WM
T2starb = 50e-3;

% Correct M0b for T2* and the partition coefficient
M0b = M0(jj) / exp(-TE/T2starref) / PC * exp(-TE/T2starb)

% Save
save('M0b.txt','M0b','-ascii');

%% Read in input data
perf = load('perf.mat'); perffnames = fieldnames(perf);
eval(['perf = perf.' perffnames{1} ';'])

%% Phase correct
perf_pc = real(perf.*exp(-1i*angle(calib_bc(:,:,:,end))));

% Check for single inversion BGS
if strcmpi(BGSMode,'SI')
% NB. We need to add a minus sign here because the single inversion pulse
% used means that the label and control conditions are switched, so the
% standard control - tag subtraction leads to a negative longitudinal
% magnetisation difference.
    perf_pc = -perf_pc;
end

% Find slices where there is reasonable perfusion signal
MaskSlcIdx = find(squeeze(sum(sum(abs(perf(:,:,:,end)) > max(abs(perf(:)))*0.1))));
slcs = round(linspace(min(MaskSlcIdx),max(MaskSlcIdx),3));

% Display pre and post phase correction to make sure there is no
% significant signal loss
DispIm([abs(perf(:,:,slcs,:)) perf_pc(:,:,slcs,:)]); title 'Absolute perfusion data (left); Phase corrected (right)'

%% Correct for coil sensitivity, M0 and PCASL inversion efficiency
PCASLInvEff = 0.85;
perf_calib = perf_pc ./ bias / M0b / PCASLInvEff; % Leave in units of s^1 for now

% Save out the calibrated perfusion data
SaveCAPRIAToNifti(perf_calib,'asldata_calib',[3.1 3.1 3.1 TempResn],[],[PerfDir '/' ReconFName '.nii.gz']);

%% Return to the original directory
cd(CurDir);