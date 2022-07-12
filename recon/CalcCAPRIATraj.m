%   Calculate the k-space trajectory associated with a specific time window
%   within a CAPRIA readout period, including the relevant raw data indices.
%
%   Tom Okell (tokell@fmrib.ox.ac.uk)
%   June 2022
%
%   Usage:
%        [LinIdx, PhsIdx, ColIdx, Phi, Theta, kspace, G, w, scale, PrepIdx] = ...
%                 CalcCAPRIATraj(Npreps, Nsegs, LinNos, EncIdx, SeqGR, ...
%                                Seq2D, MtxSz, FullMtxSz, S, UseSongOrdering, ...
%                                UseVarSpokes, UsePairedEncs, NumberOfEncCycles, ...
%                                ProjectionsPerFrame, KBSize, IsHalfSpoke)
%
%   Outputs:
%       LinIdx  =    Line number indices required to extract raw k-space data
%       PhsIdx  =    Cardiac phase indices (used for the time dimension) 
%                    required to extract raw k-space data
%       ColIdx  =    Column (readout point) indices required to extract raw
%                    k-space data
%       Phi     =    Azimuthal angles of each acquired spoke (rads)
%       Theta   =    Polar angles of each acquired spoke (rads)
%       kspace  =    kspace sampling points, size samples x kx/ky/<kz>
%                    (rads/voxel)
%       G       =    Gnufft operator corresponding to this trajectory (only
%                    generated if nargout > 6)
%       w       =    density compensation weights calculated using the Pipe
%                    method (only calculated if nargin > 7)
%       scale   =    Gnufft scaling factor (only calculated if nargout > 8)
%       PrepIdx =    Index of the ASL preparation number associated with
%                    each spoke (not required for raw data access, but
%                    useful for debugging)
%
%   Required Inputs:
%       Npreps              =  Number of ASL preparations acquired or desired 
%                              for the reconstruction. If provided as a two
%                              element vector [Npreps NprepsOffset], the
%                              offset allows some initial data to be
%                              skipped
%       Nsegs               =  Number of "segments" in the raw data,
%                              equivalent to the number of spokes acquired
%                              within the nominal temporal resolution,
%                              t_max
%       LinNos              =  Spoke numbers, relative to the start of the
%                              readout period, that should be included in
%                              the reconstruction of this frame
%       EncIdx              =  Encoding index (only required if using
%                              variable spokes across encoding cycles,
%                              otherwise set to 1)
%       SeqGR               =  Is the sequence golden ratio?
%       Seq2D               =  Is the sequence 2D (rather than 3D)?
%       MtxSz               =  The number of k-space points to be included
%                              in the reconstruction (including readout
%                              oversampling)
%       FullMtxSz           =  Number of k-space samples in a fully sampled
%                              spoke (i.e. as if partial Fourier/asymmetric
%                              echo was not used), including readout
%                              oversampling
%       S                   =  The actual size of the raw data (from
%                              twix_obj.image.dataSize), which accounts for
%                              e.g. partial Fourier/asymmetric echo
%       UseSongOrdering     =  Was Song (MRM 2014) ordering used?
%       UseVarSpokes        =  Were the acquired spokes varied across
%                              label/control/encoding cycles?
%       UsePairedEncs       =  For VarSpokes, was the same trajectory used
%                              for pairs of encoding cycles?
%       NumberOfEncCycles   =  Total number of label/control/encoding
%                              cycles (only relevant for variable spokes)
%       ProjectionsPerFrame =  Total number of spokes acquired in any t_max
%                              window across all ASL preparations in the
%                              full raw data
%
%   Optional Inputs - parameters passed as empty will take the default values:
%       KBSize              =  If the G operator is returned, this sets the
%                              size of the kaiser-bessel kernel used within
%                              the NUFFT
%       IsHalfSpoke         =  For radial-out (half spoke) data (limited
%                              testing - use with care)
%       
%
function [LinIdx, PhsIdx, ColIdx, Phi, Theta, kspace, G, w, scale, PrepIdx] = ...
                        CalcCAPRIATraj(Npreps, Nsegs, LinNos, EncIdx, SeqGR, ...
                                       Seq2D, MtxSz, FullMtxSz, S, UseSongOrdering, ...
                                       UseVarSpokes, UsePairedEncs, NumberOfEncCycles, ProjectionsPerFrame, ...
                                       KBSize, IsHalfSpoke)

    % Optional arguments
    if nargin < 16 || isempty(KBSize); KBSize = 6; end
    if nargin < 17 || isempty(IsHalfSpoke); IsHalfSpoke = false; end
    
    % Check for non-GR 3D trajectory - not yet implemented
    if ~Seq2D && ~SeqGR
        error('Non-golden ratio 3D trajectory not yet implemented!')
    end
    
    % Check if a prep offset is defined via Npreps
    if length(Npreps) > 1
        NprepsOffset = Npreps(2);
        Npreps = Npreps(1);
    else
        NprepsOffset = 0;
    end
    
    % Define some useful numbers for later
    CGR = -(1-sqrt(5))/2; % Conjugate golden ratio, Chan, MRM 2009
    PhiInc = pi/ProjectionsPerFrame; % For linear increments

    % Initialise output indices
    LinIdx = []; PhsIdx = []; PrepIdx = [];
    
    % Loop through the ASL preps to be included, adding appropriate indices
    NprepStart = NprepsOffset + 1;
    for jj = NprepStart:(NprepStart+Npreps-1)
        
        % Line number in the raw data at the start of the readout for this
        % ASL prep
        StartLinNo = (jj-1)*Nsegs + 1;
        
        % The sequence loops across "lines", then cardiac phases within a
        % readout after each ASL prep, then continues with the line
        % increments within a cardiac phase after the next ASL prep e.g.
        % for Nsegs = 3, Npreps = 2, the line numbers would be:
        %             |-Phs 1-|-Phs 2-|...
        % ASL prep 1: |-1-2-3-|-1-2-3-|...
        % ASL prep 2: |-4-5-6-|-4-5-6-|...
        % So if we wanted e.g. a temporal window that included spokes 3 and 4
        % after the ASL prep, we should extract:
        % (Lin 3, Phs 1), (Lin 1, Phs 2), (Lin 6, Phs 1) and (Lin 4, Phs 2)
        % Note we use Matlab indexing here (starting from 1), whereas the
        % sequence code indexes from 0, so this must be accounted for below
        LinIdx  = [LinIdx (StartLinNo+mod(LinNos-1,Nsegs))];
        PhsIdx  = [PhsIdx (floor((LinNos-1)/Nsegs)+1)];
        PrepIdx = [PrepIdx ones(1,length(LinNos))*jj];
    end
    
    % Calculate the appropriate azimuthal and polar angles
    if ~SeqGR % Standard equally spaced radial ordering
        Phi = (LinIdx - 1) * PhiInc;
        
        % If odd number of views, the sequence doubles the angles
        if mod(length(Phi),2)~=0
            Phi = Phi *2;
        end
        
        Theta = []; % 3D non-GR not implemented - assume 2D here
        
    else % Golden angle
        
        % Calculate a golden ratio counter (m) for each spoke
        GRCounter = calculateGRCounter(LinIdx, PhsIdx, EncIdx, UseSongOrdering, UseVarSpokes, UsePairedEncs, NumberOfEncCycles, Nsegs, ProjectionsPerFrame);

        % Calculate the relevant azimuthal and polar angles
        if Seq2D
            Phi = mod( GRCounter * CGR * pi, 2*pi );
            Theta = []; % No polar angles for 2D
        else % 3D
            [Phi, Theta] = GoldenRatio3DAngles(GRCounter,true);
        end
    end
    
    % Construct the trajectory
    [kspace, ColIdx] = CalcRadialKSpaceTraj(Phi,MtxSz,Theta,S(1),FullMtxSz,IsHalfSpoke);
    
    % Change signs to get final image in the right orientation
    if Seq2D
        kspace = [-kspace(:,1) kspace(:,2)];
    else
        kspace = [-kspace(:,1) kspace(:,2) kspace(:,3)];
    end
    
    % Also set up forward operators etc. if needed
    if nargout > 6
        
        % Set up the NUFFT
        if Seq2D
            G = Gnufft({kspace, [MtxSz MtxSz], [1 1]*KBSize, 2*[MtxSz MtxSz], [MtxSz/2 MtxSz/2]});
        else % 3D
            G = Gnufft({kspace, [MtxSz MtxSz MtxSz], [1 1 1]*KBSize, 2*[MtxSz MtxSz MtxSz], [MtxSz/2 MtxSz/2 MtxSz/2]});
        end
        
    end
    
    if nargout > 7
        % Calculate the weights using the iterative Pipe method (MRM 1999)
        % TODO: initialise with values that are proportional to abs(k) to help
        % convergence with fewer iterations?
        disp('Calculating weights')
        P = G.arg.st.p;
        w = ones([size(kspace,1) 1]);
        for kk=1:20 % NB. Hardwired to 20 iterations here
            disp(['Iteration ' ns(kk)])
            tmp = P * (P' * w);
            w = w ./ real(tmp);
        end
        
        % Clear unnecessary variables to save on memory
        clear P tmp
    end
    
    if nargout > 8
        % Calculate the scaling factor
        if Seq2D
            scale = G.arg.st.sn(end/2,end/2)^(-2) / prod(G.arg.st.Kd);
        else
            scale = G.arg.st.sn(end/2,end/2,end/2)^(-2) / prod(G.arg.st.Kd);
        end
    end