% This function saves raw CAPRIA (and other) images, RawIms, to Nifti
% format, using the filename OutFileName. If complex, separate magnitude
% and phase images are saved. Scales are the voxels sizes and temporal
% resolution, CentreSliceOffset defines the physical location of the
% imaging region relative to isocentre (in mm, only deals with
% translations, not oblique imaging slices), CopyGeomFromFName defines a
% file name to copy geometry information (voxel sizes, orientation, FOV
% offset etc.) from (overrides Scales and CentreSliceOffset).
%
% Tom Okell, June 2022
%
% [OutMagFName OutPhFName] = SaveCAPRIAToNifti(RawIms, OutFileName, Scales, ...
%                                              CentreSliceOffset, CopyGeomFromFName)

function [OutMagFName, OutPhFName] = SaveCAPRIAToNifti(RawIms, OutFileName, Scales, ...
                                                       CentreSliceOffset, CopyGeomFromFName)
  
% Deal with optional arguments
if (nargin < 3) || (sum(isnan(Scales))>0); Scales = [1 1 1 1]; end
if nargin < 4; CentreSliceOffset = [0 0 0]; end
if nargin < 5; CopyGeomFromFName = []; end

% Concatenate label/control/encoding cycles in the fifth/later dimensions
disp('Concatenating cycles...')
OutIm = reshape(RawIms,size(RawIms,1),size(RawIms,2),size(RawIms,3),[]);

% Define a permutation (Matlab rows -> y, columns -> x)
PermOrder = [2 1 3 4 5];
    
% Perform the permutation
OutIm = permute(OutIm,PermOrder);
Scales = Scales(PermOrder(1:4));
    
disp(['Size of output image is now: ' ns(size(OutIm))])
disp(['Scales are now: ' ns(Scales(:)')])
    
% Flip dimensions to ensure consistent representation 
disp('Flipping y dimension...')
OutIm = flip(OutIm,2);    
disp('Flipping z dimension...')
OutIm = flip(OutIm,3);

% Check for complex data
if ~isreal(OutIm)
    
    % Save the image using Siemens phase convention
    disp('Saving magnitude and phase images...')
    
    % Determine the output file names:
    % If the output name is <Num>_description then insert mag and ph in the
    % right place
    OutName = regexprep( OutFileName, '.*/', '' );
    DirName = regexprep( OutFileName, OutName, '');
    
    if ~isempty( regexp( OutName, '^[0-9]*_' ) )
        OutMagFName = [DirName regexprep( OutName, '_', '_mag_', 'once' )];
        OutPhFName  = [DirName regexprep( OutName, '_', '_ph_' , 'once' )];
    else % Otherwise, put mag and ph at the end of the file name
        OutMagFName = [OutFileName '_mag'];
        OutPhFName = [OutFileName '_ph'];
    end
    
    disp(['Output file names are: ' OutMagFName ' and ' OutPhFName ]);
    
    Save_Mag_Ph(OutIm, OutMagFName, OutPhFName, Scales)
    
    % Apply the correct transformation matrices
    if isempty(CopyGeomFromFName)
        tosetniftiorientinfo(OutMagFName,Scales,size(OutIm),CentreSliceOffset);
        tosetniftiorientinfo(OutPhFName,Scales,size(OutIm),CentreSliceOffset);
    else
        cmd = ['fslcpgeom ' CopyGeomFromFName ' ' OutMagFName ' -d'];
        disp(cmd); system(cmd);
        cmd = ['fslcpgeom ' CopyGeomFromFName ' ' OutPhFName  ' -d'];
        disp(cmd); system(cmd);
    end
    
else
    disp('Data not complex, saving real values...')
    save_avw( OutIm, OutFileName, 'f', Scales );
    
    % Apply the correct transformation matrices
    if isempty(CopyGeomFromFName)
        tosetniftiorientinfo(OutFileName,Scales,size(OutIm),CentreSliceOffset);
    else
        cmd = ['fslcpgeom ' CopyGeomFromFName ' ' OutFileName ' -d'];
        disp(cmd); system(cmd);
    end
end