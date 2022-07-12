% Calculate a golden ratio counter for a a given line, phase, encoding
% index for various sequence parameters, mimicking the CAPRIA sequence
% behaviour
function GRCounter = calculateGRCounter(LinIdx, PhsIdx, EncIdx, UseSongOrdering, UseVarSpokes, UsePairedEncs, NumberOfEncCycles, Segments, ProjectionsPerFrame)

% Convert Matlab to C++ array indexing
LinIdx = LinIdx - 1;
PhsIdx = PhsIdx - 1;
EncIdx = EncIdx - 1;

% Convert nominal encoding indices to actual encoding indices
if UseVarSpokes
    if UsePairedEncs        
        EncIdx = floor( EncIdx / 2 );        
        % NB. For unpaired spokes, the unique encoding number is the same as the actual encoding number
        % so no need to make any modifications
    end
    
else % No variation across encoding cycles, so the unique encoding number is always zero
    EncIdx(:) = 0;    
end

% Determine the number of encoding cycles with unique spokes
if UseVarSpokes
    if UsePairedEncs
        % For paired cycles, we need an even number of encodings
        if mod(NumberOfEncCycles,2) ~= 0
            error(">>> Error - paired encodings cannot be used with an odd number of encoding cycles");
        end
        
        NumberOfUniqueEncCycles = NumberOfEncCycles / 2;
        
    else % Variable across spokes, unpaired, so all encoding cycles have unique spokes        
        NumberOfUniqueEncCycles = NumberOfEncCycles;        
    end
    
else % No variation across encoding cycles
    NumberOfUniqueEncCycles = 1;
end


% Determine the number of ASL preps for the entire scan (per VEPCASL cycle)
NoVEPCASLPrepsWithUniqueSpokes = ceil( ProjectionsPerFrame / Segments ) * NumberOfUniqueEncCycles;

% Work out which ASL prep we are on (for variable spokes, NumberOfUniqueEncCycles will be more than one and
% lEnc will go above zero so we also include different encodings)
CurrentVEPCASLPrep = floor( (LinIdx / Segments) ) * NumberOfUniqueEncCycles + EncIdx;

% Find the line number (as set in the MDH) corresponding to the first line acquired in phase 0 for this ASL prep
StartLinNo = floor( LinIdx / Segments ) * Segments;

% How far is the current line from the first line acquired in this readout block?
DistFromROStart = LinIdx + PhsIdx * Segments - StartLinNo;

% Determine the golden ratio counter, which will depend on the looping required
if (UseSongOrdering) % Order spokes as per Song et al MRM 2014, which increase down repeats then across time
    GRCounter = CurrentVEPCASLPrep + DistFromROStart * NoVEPCASLPrepsWithUniqueSpokes;
    
else % Not Song ordering: default to tmax method (Okell MRM 2019)
    % T.O. If the "TR" (i.e. temporal resolution) is set to the maximum desired temporal resolution
    % (t_max) then for each new ASL prep, the line number will continue counting from this maximum value.
    % The number of spokes within a given interval (and thus the total scan time) can be tuned using
    % the "radial views" option.
    GRCounter = DistFromROStart + CurrentVEPCASLPrep * Segments;
    
    % Old version for same spokes at each encoding: the equations above give the same result when m_lNumberOfUniqueEncCycles = 1, which
    % is the case when variable spokes is switched off.
    %	GRCounter = LinIdx + PhsIdx * Segments;
end

