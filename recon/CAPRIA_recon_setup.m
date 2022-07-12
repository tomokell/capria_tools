% Add the paths required for CAPRIA reconstructions
%
% Tom Okell, June 2022

disp('>> CAPRIA recon setup starting...')

%% Add the directory containing this script and subdirectories
% If we're running the full script, mfilename should work:
filePath = mfilename('fullpath');

% This could have failed if running separate sections of the script in the Matlab editor, so try an alternative approach
% NB. Running section by section seems to result in mfilename returning a
% temporary file in /private/var, which doesn't contain the other mfiles
% needed, so in that case also try this alternative approach
if isempty(filePath) || ~isempty(regexp(filePath,'/private/var/', 'once'))
    filePath = matlab.desktop.editor.getActiveFilename;
end
if isempty(filePath) % Error out if we can't find the location of this script
    error('Could not locate the current script')
end
disp(['Current script is: ' filePath]);

% Find the containing directory
filePath = fileparts(filePath);
disp(['Adding ' filePath ' and subdirectories to the Matlab path...']);

addpath(genpath(filePath)); % Adds all the folders and subfolders within the containing folder to the matlab path

% Check if there is a utils directory at the same level in the directory
% structure
utilsdir = [filePath '/../utils'];
if exist(utilsdir,'dir')
    disp(['Also adding ' utilsdir]);
    addpath(utilsdir)
end

%% Add the FSL matlab directory (for reading/writing Nifti files)
fsldir = getenv('FSLDIR'); % Find the FSL directory from the FSLDIR environment variable
% If the environment variable isn't set properly, you can specify the FSL
% directory explicitly like this:
% fsldir = '/usr/local/fsl';

if isempty(fsldir) % If it's not set, see if the default location exists
    if exist('/usr/local/fsl','dir')
        fsldir = '/usr/local/fsl';
    else
        error(['FSLDIR environment variable not specified and default location not found: ' ...
            'please ensure FSL is installed and FSLDIR is set in the terminal prior to calling Matlab ' ...
            'or edit the calling Matlab script to specify the FSL directory explicitly']);
    end
end
disp(['Found FSL directory: ' fsldir]);
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
disp(['Adding ' fsldirmpath ' to the Matlab path...']);
addpath(fsldirmpath);


%% Add the IRT tools (including the NUFFT)
irtsetup = which('irt/setup');
if isempty(irtsetup)
    error('Cannot find IRT setup script - please ensure IRT is installed and on the Matlab path')
end

% Set the irtdir variable, which is used by irtsetup
irtdir = fileparts(irtsetup);
disp(['Found IRT directory: ' irtdir])

% Running
disp(['Running IRT setup file: ' irtsetup]);
run(irtsetup)

%% Done!
CAPRIA_recon_setup_complete = true;
disp('>> CAPRIA recon setup completed successfully!')

