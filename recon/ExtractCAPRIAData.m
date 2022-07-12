%   Extract combined angiography and perfusion using radial imaging and ASL 
%   (CAPRIA) data from a Siemens raw data (meas.dat) file. Requires
%   mapVBVD.m.
%
%   Tom Okell (tokell@fmrib.ox.ac.uk)
%   June 2022
%
%   Usage:
%        [kdata, kspace, S, CentreSlcOffset, Resn] = ...
%           ExtractCAPRIAData(MeasFileID, Seq2D, TempResn, FracDataToUse, ...
%           subtract, SpatialMtxSz, SeqType, PhaseCorr, FrameNos, ...
%           ConcatROPointsAndSpokes, IsHalfSpoke)
%
%   Outputs:
%       kdata            =    k-space data 
%                             dims: [samples <spokes> time coils
%                                    label/control])
%       kspace           =    k-space sample positions (in rads/voxel)
%                             dims: [samples <spokes> time kx/ky/<kz>
%       S                =    size of the raw data stored in meas.dat
%       CentreSlcOffset  =    3D vector containing the physical
%                             displacement of the centre of the FOV
%                             from isocentre (mm)
%       Resn             =    4D vector containing the voxel dimensions and
%                             temporal resolution (mm and s)
%
%   Required Inputs:
%       MeasFileID    =    Numeric measurement ID of the meas.dat file
%       Seq2D         =    Is the sequence 2D? Leave empty to determine
%                          this from the raw data file
%       TempResn      =    Temporal resolution (ms) to reconstruct at
%
%   Optional Inputs - parameters passed as empty will take the default values:
%       FracDataToUse =    Fraction of data to use in the reconstruction.
%                          Specify as a single value (e.g. 0.5 would use
%                          the first half of the raw data) or a two
%                          component vector [Frac Offset], to specify the
%                          fraction of data to use and an offset from the
%                          beginning of the file (e.g. [0.5 0.5] would use
%                          the second half the data - half the data, 
%                          starting from halfway through).
%       subtract      =    subtract label and control data in k-space (if
%                          false, both label and control data are returned)
%       SpatialMtxSz  =    nominal matrix size to reconstruct at
%                          (determines reconstructed voxels size via
%                          nominal FOV (from raw data) / SpatialMtxSz
%       SeqType       =    CAPRIA sequence type: default is to read
%                          from the raw data header, which works for
%                          all current CAPRIA sequences.
%       PhaseCorr     =    perform phase correction step to align
%                          k-space data acquired in label/control
%                          conditions
%       FrameNos      =    Frame numbers to reconstruct (leave empty for
%                          all)
%       ConcatROPointsAndSpokes =    Reshape the output data so that readout points
%                          and spokes are all combined into the first
%                          dimension. If true, output dims of kdata are: 
%                          [kspace_points time coils label/control/encodings]
%                          If false, output dims are:
%                          [kspace_points spokes time coils label/control/encodings]
%       IsHalfSpoke   =    For radial-out (half spoke) data (limited
%                          testing - use with care)
%       
%
function [kdata, kspace, S, CentreSlcOffset, Resn] = ExtractCAPRIAData(MeasFileID, Seq2D, TempResn, FracDataToUse, subtract, SpatialMtxSz, SeqType, PhaseCorr, FrameNos, ConcatROPointsAndSpokes, IsHalfSpoke)

if nargin < 3  || isempty(TempResn);        TempResn = [];          end
if nargin < 4  || isempty(FracDataToUse);   FracDataToUse = 1;      end
if nargin < 5  || isempty(subtract);        subtract = true;        end
if nargin < 6  || isempty(SpatialMtxSz);    SpatialMtxSz = [];      end
if nargin < 7  || isempty(SeqType);         SeqType = [];           end
if nargin < 8  || isempty(PhaseCorr);       PhaseCorr = true;       end
if nargin < 9  || isempty(FrameNos);        FrameNos = [];          end
if nargin < 10 || isempty(ConcatROPointsAndSpokes);   ConcatROPointsAndSpokes = true;   end
if nargin < 11 || isempty(IsHalfSpoke);     IsHalfSpoke = false;    end

% Read the raw data headers
twix_obj = mapVBVD(MeasFileID,'ignoreSeg',true,'removeOS',false);

% If there are pre-scan adjustment data, discard here
if numel(twix_obj) > 1
    twix_obj = twix_obj{end};
end

% Extract the sequence name
SeqName = twix_obj.hdr.Config.SequenceFileName;
disp(['Sequence used was: ' SeqName])

% Determine the sequence parameters based on the sequence version. We want
% to determine:
% Seq2D         = Is the sequence 2D (i.e. not 3D)?
% SeqGR         = Was the golden ratio trajectory used?
% SeqSong       = Was Song et al. (MRM 2014) sequence looping used?
% SeqVarSpokes  = Were the spokes variable across label/control/encoding conditions
% SeqPairedEncs = Were the spokes the same for pairs of encoding cycles
%                 (used for vessel encoded data only)
switch SeqName
    case {'%CustomerSeq%\to_CAPIASL_CV_nce_angio','%CustomerSeq%\to_CV_VEPCASL'} % Original version
        
        % Determine the data type: 1 = linear, 2 = 2D GR, 3 = 3D GR
        if isempty(SeqType)
            SeqType = twix_obj.hdr.MeasYaps.sWiPMemBlock.alFree{12+1};
        else
            % SeqType is provided - useful for old data where the above parameter isn't
            % defined in the same way.
        end
        
        switch SeqType
            case 1
                Seq2D = true;   
                SeqGR = false;
                SeqSong = false;
                
            case 2
                Seq2D = true;
                SeqGR = true;
                SeqSong = false;
                
            case 4
                Seq2D = false;
                SeqGR = true;
                SeqSong = false;
                
            case 24
                Seq2D = true;
                SeqGR = true;
                SeqSong = true;
                
            otherwise
                error(['SeqType = ' ns(twix_obj.hdr.MeasYaps.sWiPMemBlock.alFree{12+1}) ': not recognised!']);
        end

        % For this older sequence, variable spokes and paired encodings were not used
        SeqVarSpokes = false; SeqPairedEncs = false;
        
    case {'%CustomerSeq%\to_CV_VEPCASL_v0p2','%CustomerSeq%\to_CV_VEPCASL_v0p3'}
        
        % WIP mem block position, as defined in the sequence code, to
        % Matlab index conversion:
        % +1 for C++ vs. Matlab indexing
        % +2 for reserved WIP mem block positioned used elsewhere in the sequence
        % +WIPMemBlock starting position (7 for trajectory boolean parameters)
        % +Position within the boolean array
        % NB. False yields an empty array, true yields 1
        Seq2D = isempty(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{1 + 2 + 7 + 0});
        if Seq2D disp('Sequence is 2D'); else disp('Sequence is 3D'); end
        
        SeqGR = ~isempty(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{1 + 2 + 7 + 1});
        if SeqGR disp('Sequence uses GR'); else disp('Sequence does not use GR'); end
        
        SeqSong = ~isempty(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{1 + 2 + 7 + 2});
        if SeqSong disp('Sequence uses Song ordering'); else disp('Sequence does not use Song ordering'); end
        
        SeqVarSpokes = ~isempty(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{1 + 2 + 7 + 3});
        if SeqVarSpokes disp('Sequence uses variable spokes'); else disp('Sequence does not use variable spokes'); end
        
        SeqPairedEncs = ~isempty(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{1 + 2 + 7 + 4});
        if SeqPairedEncs disp('Sequence uses paired encodings'); else disp('Sequence does not use paired encodings'); end

case {'%CustomerSeq%\to_CV_VEPCASL_v0p4'}
        
        % WIP mem block positions as above, but now WIPMemBlock starting position 
        % is 8 for trajectory boolean parameters
        Seq2D = isempty(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{1 + 2 + 8 + 0});
        if Seq2D disp('Sequence is 2D'); else disp('Sequence is 3D'); end
        
        SeqGR = ~isempty(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{1 + 2 + 8 + 1});
        if SeqGR disp('Sequence uses GR'); else disp('Sequence does not use GR'); end
        
        SeqSong = ~isempty(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{1 + 2 + 8 + 2});
        if SeqSong disp('Sequence uses Song ordering'); else disp('Sequence does not use Song ordering'); end
        
        SeqVarSpokes = ~isempty(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{1 + 2 + 8 + 3});
        if SeqVarSpokes disp('Sequence uses variable spokes'); else disp('Sequence does not use variable spokes'); end
        
        SeqPairedEncs = ~isempty(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{1 + 2 + 8 + 4});
        if SeqPairedEncs disp('Sequence uses paired encodings'); else disp('Sequence does not use paired encodings'); end

case {'%CustomerSeq%\to_CV_VEPCASL_v0p6','%CustomerSeq%\to_CV_VEPCASL_v0p5','%CustomerSeq%\qijia_CV_VEPCASL_v0p4','%CustomerSeq%\qijia_CV_VEPCASL'}
        
        % WIP mem block positions as above, but now WIPMemBlock starting position 
        % is 9 for trajectory boolean parameters
        Seq2D = isempty(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{1 + 2 + 9 + 0});
        if Seq2D disp('Sequence is 2D'); else disp('Sequence is 3D'); end
        
        SeqGR = ~isempty(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{1 + 2 + 9 + 1});
        if SeqGR disp('Sequence uses GR'); else disp('Sequence does not use GR'); end
        
        SeqSong = ~isempty(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{1 + 2 + 9 + 2});
        if SeqSong disp('Sequence uses Song ordering'); else disp('Sequence does not use Song ordering'); end
        
        SeqVarSpokes = ~isempty(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{1 + 2 + 9 + 3});
        if SeqVarSpokes disp('Sequence uses variable spokes'); else disp('Sequence does not use variable spokes'); end
        
        SeqPairedEncs = ~isempty(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{1 + 2 + 9 + 4});
        if SeqPairedEncs disp('Sequence uses paired encodings'); else disp('Sequence does not use paired encodings'); end
                        
    otherwise
        disp(['Unknown sequence version: ' SeqName ', continuing with provided parameters'])
        
end

% Variable spokes not yet supported, so throw an error here
if SeqVarSpokes
    error('Variable spokes not yet supported!')
end

% Extract useful info from the header
MaxTempResn = twix_obj.hdr.Config.TR/1000; % t_max, in ms
Nsegs = twix_obj.hdr.Config.NSeg;          % Number of excitations in t_max (M)
TotNSpokesInMaxTempResn = twix_obj.hdr.Config.RawLin; % total number of spokes acquired in any t_max window across all ASL preparations
NPhases = twix_obj.hdr.Config.NPhs; % Number of phases (frames, if temporal resolution is set to t_max
FullMtxSz = twix_obj.hdr.Config.RawCol; % Number of k-space samples in a full spoke, including readout oversampling

% If not provided, set temporal resolution of the recon to t_max
if isempty(TempResn)
    disp(['Temporal resolution not provided, setting to t_max: ' ns(MaxTempResn) ' ms'])
    TempResn = MaxTempResn;
end

% If not provided, set the spatial matrix size to be the maximum,
% accounting for readout oversampling
if isempty(SpatialMtxSz)    
    SpatialMtxSz = FullMtxSz/2;
    disp(['Spatial resolution not provided, setting to maximum: ' ns(SpatialMtxSz)])
end

% For non-golden ratio reconstructions, the temporal resolution is fixed
if (SeqGR == false) && (TempResn ~= MaxTempResn)
    warning(['Specified temporal resolution (' ns(TempResn) ') does not match the value set in the protocol: ' ns(MaxTempResn)])
    warning(['Reverting to protocol temporal resolution: ' ns(MaxTempResn)])
    TempResn = MaxTempResn;
end

% Determine the data points to use for reconstruction
% Find the index of the k=0 point
kcentreIdx = ceil((FullMtxSz+1)/2);

% Round the spatial matrix size to the nearest even number
if mod(round(SpatialMtxSz),2) ~=0    
    SpatialMtxSz = round(SpatialMtxSz)+1;
    disp(['Rounding spatial matrix size up to the nearest even number: ' (SpatialMtxSz)])
end

% Specify the number of k-space samples to be used from each spoke
MtxSz = SpatialMtxSz * 2; % Account for RO oversampling

% Define an index for extracting data from each spoke
ColIdx = (kcentreIdx-MtxSz/2):(kcentreIdx+MtxSz/2-1);

% Calculate useful parameters
TR = MaxTempResn / Nsegs; % Time between excitations, ms
Npreps = TotNSpokesInMaxTempResn / Nsegs; % Number of ASL preps per label/control/encoding cycle
disp([ns(Npreps) ' ASL preparations were used per encoding cycle'])
Nspokes = round(TempResn/TR); % Number of spokes acquired within the temporal resolution after each ASL prep
NewTempResn = TR*Nspokes; % Round temporal resolution to an integer number of TRs
disp(['Rounding temporal resolution to ' ns(NewTempResn)]);
NSpokesPerPrep = Nsegs * NPhases; % Number of spokes acquired after each ASL prep
T = NPhases * MaxTempResn; % Total readout time
disp(['Total readout time was: ' ns(T) ' ms'])
NFrames = floor(NSpokesPerPrep/Nspokes); % Number of frames to reconstruct
disp(['This yields ' ns(NFrames) ' frames']);
VoxSz = twix_obj.hdr.Config.ReadFoV / MtxSz * 2; % Voxel size, accounting for RO oversampling

% Determine the slice thickness for 2D/3D data separately
if Seq2D
    SlcThk = twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.dThickness;
else
    SlcThk = VoxSz;
end
S = twix_obj.image.dataSize; % Grab the input data size

% Save out the voxel size (in mm) and temporal resolution (in s)
Resn = [VoxSz VoxSz SlcThk NewTempResn/1000];

% Save out the central slice offset
if isfield(twix_obj.hdr.Phoenix.sSliceArray.asSlice{1},'sPosition')
    SlicePos = twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.sPosition;
else
    disp('Could not find central slice location: assuming isocentre')
    SlicePos = [];
end

if isfield(SlicePos,'dSag'); Sag = SlicePos.dSag; else; Sag = 0; end
if isfield(SlicePos,'dCor'); Cor = SlicePos.dCor; else; Cor = 0; end
if isfield(SlicePos,'dTra'); Tra = SlicePos.dTra; else; Tra = 0; end

CentreSlcOffset = [Sag Cor Tra];

% Cut the number of preps to use for recon if necessary
if length(FracDataToUse) > 1 % If we want to use less than 100% of the data with an offset, deal with that here
    FracDataOffset = FracDataToUse(2);
    FracDataToUse = FracDataToUse(1);
    
    % Check this makes sense (i.e. we don't try to read past the end of the
    % data
    if FracDataToUse + FracDataOffset > 1
        error(['FracDataToUse(' ns(FracDataToUse) ') + FracDataOffset (' ns(FracDataOffset) ') is greater than 1!!'])
    end
else
    FracDataOffset = 0;
end

% If we want to offset, calculate this first
NprepsOffset = round(Npreps * FracDataOffset);

% Now reset Npreps to be reduced if needed
if FracDataToUse ~=1
    disp(['Only using ' ns(round(Npreps * FracDataToUse)) ' of ' ns(Npreps) ' ASL preparations']);  
    disp(['Using an offset of ' ns(NprepsOffset) ' preps']);
    Npreps = round(Npreps * FracDataToUse);
end

% Determine the acceleration (undersampling) factor, R
if Seq2D % 2D radial
    R = SpatialMtxSz * pi/2 / (Nspokes*Npreps);    
else % 3D radial    
    R = (SpatialMtxSz)^2 * pi/2 / (Nspokes*Npreps);    
end
disp(['Undersampling factor = ' ns(R)])
    
% If no frame numbers are specified, reconstruct all
if isempty(FrameNos)
    FrameNos = 1:NFrames;
end
    
% Assume the number of encoding cycles is the number of measurements * averages 
% for now since both of these loops increment the encoding cycle number.
% Note that for standard CAPRIA with PCASL label and control, this should
% be 2, but we try to keep this general for vessel or time-encoded data
% where more encoding cycles are necessary
NumberOfEncCycles = twix_obj.image.dataSize(6) * twix_obj.image.dataSize(9);
disp(['Averages in raw data: ' ns(twix_obj.image.dataSize(6)) ...
      ', Measurements in raw data: ' ns(twix_obj.image.dataSize(9))])
disp(['Assuming we have ' ns(NumberOfEncCycles) ' encoding cycles'])

% Check if we want subtraction, we need an even number of encoding cycles
if subtract && (mod(NumberOfEncCycles,2) ~= 0)
    error('Cannot use subtraction with an odd number of encoding cycles!')
end

% Initialise the output
kspace=[];

% Grab the raw data with the appropriate indices
disp('Reading in raw data...')

% Loop through the output frames
for ii = 1:length(FrameNos)
    
    disp(['Getting kspace trajectory and data for frame ' ns(ii) ' of ' ns(length(FrameNos))])
    
    NFrame = FrameNos(ii);
    
    % Identify the actual spoke numbers relative to the start of the
    % readout that we want to grab here
    LinNos = ((NFrame-1)*Nspokes+1):(NFrame*Nspokes);
    
    % Empty the index of relevant line and phase numbers
    LinIdx = []; PhsIdx = [];
    
    % Assume we are not using variable spokes for now, so all trajectories are the same for each encoding
    EncIdx = 1; 
    
    % Find the raw data indices we want and the corresponding k-space
    % sample locations
    [LinIdx, PhsIdx, ColIdx, Phi, ~, kspace(:,ii,:)] = CalcCAPRIATraj([Npreps NprepsOffset], Nsegs, LinNos, EncIdx, SeqGR, Seq2D, MtxSz, FullMtxSz, S, ...
                                                                      SeqSong, SeqVarSpokes, SeqPairedEncs, NumberOfEncCycles, ...
                                                                      TotNSpokesInMaxTempResn,[],IsHalfSpoke);
    
    % Loop through the radial spokes and extract the relevant data from the
    % twix file
    for kk = 1:length(Phi)
        % Dimensions of kdata: RO_cols  Coils  spokes  averages echoes measurements time
        kdata(:,:,kk,:,:,:,ii) = twix_obj.image(ColIdx,:,LinIdx(kk),:,:,:,PhsIdx(kk),:,:,:,:,:,:,:,:,:);
    end
end

% Reorganise the output data
Nt  =   size(kdata,7); % Number of frames being reconstructed
kdata = reshape(permute(kdata,[1,3,7,2,4,5,6]),length(ColIdx), [], size(kdata,2), NumberOfEncCycles);
% Dimensions of kdata now: RO_cols  spokes/time Coils averages/measurements

% Perform phase correction if requested
if PhaseCorr
    
    % Loop through encoding cycles
    for jj = 1:size(kdata,4)
        
        % Loop through all spokes and time points
        for ii = 1:size(kdata,2)
            
            % Calculate the phase correction angle required (Schauman, MRM
            % 2020) relative to the first encoding
            p=angle(mean(reshape(kdata(:,ii,:,jj).*conj(kdata(:,ii,:,1)),[],1)));
            
            % Apply it to the jth encoding cycle data to align it to the
            % first encoding
            kdata(:,ii,:,jj)=kdata(:,ii,:,jj).*exp(-1j*p);
        end
    end
end

% Subtract pairs of encoding cycles if requested
if subtract
    kdata  =   kdata(:,:,:,2:2:end)-kdata(:,:,:,1:2:end);
    EncDimSz = NumberOfEncCycles/2; % This dimension is now smaller
    EncTxt = ' Subtracted_Encoding_Cycles';
else % Do not subtract label/control/encoded data
    disp('No Subtraction');
    EncDimSz = NumberOfEncCycles;
    EncTxt = ' Encoding_cycles';
end

% Reshape the output
if ConcatROPointsAndSpokes % Concatenate readout points and spokes
    kdata   =   reshape(kdata, [],Nt,size(kdata ,3),EncDimSz);
    kspace  =   reshape(kspace,[],Nt,size(kspace,3));
    ROSpokesTxt = 'kspace_samples';
    
else % Keep spoke dimension easily accessible e.g. for debugging
    kdata   =   reshape(kdata, length(ColIdx),[],Nt,size(kdata ,3),EncDimSz);
    kspace  =   reshape(kspace,length(ColIdx),[],Nt,size(kspace,3));
    ROSpokesTxt = 'kspace_samples spokes';
    
end

disp(['Dimensions of kdata:  (' regexprep(ns(size(kdata)),'\s+',',')  '): ' ROSpokesTxt ' time Coils' EncTxt]);
disp(['Dimensions of kspace: (' regexprep(ns(size(kspace)),'\s+',',') '): ' ROSpokesTxt ' time kx/ky/<kz>']);
