% Fit CAPRIA angiographic data using a custom version of fabber, which must
% be installed and compiled, as per the main fabber instructions, including
% the modifications here (the pcasldisp branch):
% https://github.com/tomokell/fabber_models_asl/tree/pcasldisp
% We have only implemented the quadratic flip angle modulation in fabber,
% so that is assumed here also.
%
% Tom Okell, July 2025
%
%   Usage: FitCAPRIAAngioData(MatFileFName,OutDir,tau,FAParams,t0,TR,Nsegs,Nphases)
%
%   Required Inputs:
%       MatFileFName    =   Name of the CAPRIA reconstruction .mat
%                           file
%       OutDir          =   Output directory
%
%   Optional inputs:
%       tau             =   Labelling duration (s)
%       FAParams      =     Flip angle parameters (see CalcCAPRIAFAs.m)
%       t0            =     Time at which imaging starts relative to the 
%                           start of PCASL labelling (s)
%       TR            =     Time between excitation pulses in the readout
%                           (s)
%       Nsegs         =     Number of excitations within the temporal
%                           resolution (i.e. TempResn / TR)
%       Nphases       =     Number of frames (i.e. total readout time /
%                           (Nsegs*TR) )
%       T1b           =     T1 of blood (s)
%       MaskThr       =     Threshold used to generate the initial mask
%                           (after rescaling phase-corrected data by 1e10)
%       ClustersToKeep=     Number of clusters to keep after initial
%                           thresholding and clustering
%       FabberFName   =     Filename of the custom FABBER executable

function FitCAPRIAAngioData(MatFileFName,OutDir,tau,FAParams,t0,TR,Nsegs,Nphases,T1b,MaskThr,ClustersToKeep,FabberFName)

if nargin < 3; tau = 1.8;  end
if nargin < 4; FAParams = [2 9];  end
if nargin < 5; t0 = tau+2e-3+10e-3+11e-3+10e-3; end
if nargin < 6; TR = 9.1e-3;  end
if nargin < 7; Nsegs = 36/2;  end
if nargin < 8; Nphases = 6*2; end
if nargin < 9; T1b = 1.65; end
if nargin < 10; MaskThr = 35; end
if nargin < 11; ClustersToKeep = 2; end
if nargin < 12; FabberFName = '~/Documents/C++/fabber_models_asl/fabber_asl'; end

% Read in the data and extract the relevant image data from the struct
img = load(MatFileFName); imgfnames = fieldnames(img);
eval(['img = img.' imgfnames{1} ';'])

% Assume there is a Nifti file also with the correct header we can copy
NiftiFName = [regexprep(MatFileFName,'\.mat','') '.nii.gz'];

% Phase correct the data
img_pc = PhaseCorrectDynAngioIms(img,4);

% Make the output directory
mkdir(OutDir);

% Rescale to avoid small number issues and save out, copying the header
% from the recon Nifti file
SaveCAPRIAToNifti(img_pc*1e10,[OutDir '/data_pc'],[1 1 1 1],[0 0 0],NiftiFName);

% Move to the output directory
CurDir = pwd;
cd(OutDir)

% Generate a mask
tosystem(['fslmaths data_pc -nan -Tmax -thr ' ns(MaskThr) ' -bin mask']);

% Cluster
tosystem('cluster --in=mask --thresh=0.1 -o mask_clusters');

% Find the intensity allocated to the Nth largest cluster
[~,tmp]=builtin('system','fslstats mask_clusters -R');
maxI = split(tmp,' ');
maxI = str2num(maxI{2});
thr = ns(maxI - (ClustersToKeep-1));
tosystem(['fslmaths mask_clusters -thr ' thr ' -bin mask_clusters_bin'])

%% Define the timing
t = t0:TR:(t0+(Nsegs*Nphases-1)*TR);
tAv = zeros(Nphases,1);
for ii = 1:Nphases
    Start = (ii-1)*Nsegs+1;
    End = ii*Nsegs;
    tAv(ii) = mean(t(Start:End));
end

%% Run the fitting
fabcmd = FabberFName;
fabcmd = [fabcmd ' --data=data_pc.nii.gz --mask=mask_clusters_bin.nii.gz'];
fabcmd = [fabcmd ' --model=aslrest --disp=gamma --method=vb --inferdisp']; 
fabcmd = [fabcmd ' --batart=0.5']; 
fabcmd = [fabcmd ' --noise=white --allow-bad-voxels --max-iterations=20 --convergence=trialmode --max-trials=10'];
fabcmd = [fabcmd ' --save-mean --save-mvn --save-std --save-model-fit --save-residuals'];
for jj = 1:length(tAv)
    fabcmd = [fabcmd ' --ti' ns(jj) '=' ns(tAv(jj)) ' --rpt' ns(jj) '=1']; 
end
fabcmd = [fabcmd ' --tau=' ns(tau) ' --casl --slicedt=0.0 --t1=1.3 --t1b=' ns(T1b) ' --bat=1.3 --batsd=10.0 --incbat --inferbat --incart --inferart '];
fabcmd = [fabcmd ' --capria --capriafa1=' ns(FAParams(1)) ' --capriafa2=' ns(FAParams(2)) ' --capriatr=' ns(TR)];
fabcmd1 = [fabcmd ' --output=fabber_out'];

tosystem(fabcmd1)

% Convert outputs
s = ra('fabber_out_latest/mean_disp1'); % Fabber fits log(s) to prevent s becoming negative
s = exp(s);
sp = ra('fabber_out_latest/mean_disp2'); % The second parameter is log(s*p)
sp = exp(sp);
sp(sp > 10) = 10; % This mimics what is done in fabber
p = sp ./ s;
save_avw(s,'fabber_out_latest/mean_disp_s','f',[1 1 1 1])
save_avw(p,'fabber_out_latest/mean_disp_p','f',[1 1 1 1])

% Copy the Nifti header
tosystem('fslcpgeom fabber_out_latest/mean_disp1 fabber_out_latest/mean_disp_s -d')
tosystem('fslcpgeom fabber_out_latest/mean_disp1 fabber_out_latest/mean_disp_p -d')

% Go back to the current directory
cd(CurDir)