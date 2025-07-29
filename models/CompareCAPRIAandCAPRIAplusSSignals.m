% Compare the evolution of the static tissue and blood mangetization for
% CAPRIA (pre-saturation only) vs. CAPRIA+S (pre-saturation and a single
% inversion pulse prior to the readout).
%
% Note that this code relies on a simple Bloch simulator which can be found
% here: https://github.com/tomokell/bloch_sim
%
% Tom Okell, July 2025

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

% Find and add the containing directory
filePath = fileparts(filePath);
disp(['Adding ' filePath ' to the Matlab path...']);
addpath(filePath);

% Check if there is a utils directory at the same level in the directory
% structure
utilsdir = [filePath '/../utils'];
if exist(utilsdir,'dir')
    disp(['Also adding ' utilsdir]);
    addpath(utilsdir)
end

% Check if the bloch simulator is available
if isempty(which('bloch_sim.m'))
    error('bloch_sim.m not available. Please install from here: https://github.com/tomokell/bloch_sim and add to the Matlab path to continue...')
end

%% Set up tissue parameters
Tis(1).Name = 'WM';    Tis(1).T1 = 830e-3;  Tis(1).T2 =   80e-3;
Tis(2).Name = 'GM';    Tis(2).T1 = 1330e-3; Tis(2).T2 =  100e-3;
Tis(3).Name = 'CSF';   Tis(3).T1 = 4000e-3; Tis(3).T2 = 2000e-3;
Tis(4).Name = 'Blood'; Tis(4).T1 = 1932e-3; Tis(4).T2 =  275e-3;  % @ 3T, ref: Stanisz, MRM 2005
Tis(5).Name = 'Blood Labelled'; Tis(5).T1 = 1932e-3; Tis(5).T2 =  275e-3;  
LHR = true; % Use LHR rotations

%% Generate an SPGR readout
% Set up parameters
set_SSFP_defaults;  % Some sensible defaults

% Modify default parameters
FlipAng = [2 9];  VFAType = 'Quad'; % Quadratic flip angle scheme
FlipAngStr = ['_FlipAng' num2str(FlipAng(1)) '_' ns(FlipAng(2))];
ORFreq  = 0; SeqType = 'SPGR'; TR = 9.1e-3; 
tau = 1.8; T = ceil(2/TR)*TR;  dt = TR/50;

% Save a string defining off-resonance (not used here)
if ORFreq == 0; ORFreqStr  = '';
else            ORFreqStr  = ['_ORFreq' num2str(ORFreq)];
end

SaveStr = [FlipAngStr ORFreqStr];

% Generate the time array
t = single(0:dt:T);

% Generate the off resonance frequency array
OffRes = single(ones(1,size(t,2))*ORFreq);

% Set up the position array of the spin
P = single(generate_position('const_vel', Ps, [0 0 0]', t));

% Generate the gradient and RF arrays for an SPGR sequence
[G, RF, Spoil] = make_SSFP_seq(FlipAng, RF_dur, TR, dt, T, G_amp, SeqType,[],VFAType);
G = single(G);
RF = single(RF);

%% Extend the readout so there is tau of time prior to the readout
tfull = single([-tau:dt:(0-dt) t]);
Nfull = size(tfull,2) - size(t,2);
Gfull = single([zeros(3,Nfull) G]);
RFfull = single([zeros(2,Nfull) RF]);
Pfull = single([zeros(3,Nfull) P]);
OffResfull = single(ones(1,size(tfull,2)) * ORFreq );
Spoilfull = [false(1,Nfull) Spoil];

%% Generate a 180 degree RF pulse
RF180dur = 1e-3;
RF180 = single(sinc_RF(180, 0, RF180dur, dt, 7, g, 0.1, false));

%% Create the CAPRIA+S sequence with a saturation followed by inversion prior to readout
RFsatinv = RFfull;
ROStart = find(RFfull(1,:) > 0, 1, 'first');
RFInvStart = ROStart - size(RF180,2) - round(RF180dur/dt);
RFsatinv(:,RFInvStart:(RFInvStart + size(RF180,2) -1)) = RF180;

%% Simulate sequences original (sat only) and CAPRIA+S (sat inv) sequences
for TisNo = 1:size(Tis,2)
    
    disp(['Simulating tissue ' num2str(TisNo)])

    % Static tissue feels the pre-saturation, but blood doesn't see it
    if strcmp(Tis(TisNo).Name, 'Blood')  % Start at M0 since fresh blood is arriving
        Ms = [0 0  1]';
    elseif strcmp(Tis(TisNo).Name, 'Blood Labelled')  % Start at -M0, assuming perfect inversion at the start of simulation
        Ms = [0 0 -1]';
    else  % Static tissue is pre-saturated
        Ms = [0 0  0]';
    end
        
    % Seq 1: Pre-saturation only
    Tis(TisNo).MSat = bloch_sim(Pfull, Gfull, RFfull, OffResfull, tfull, Tis(TisNo).T1, Tis(TisNo).T2, Ms, LHR, Spoilfull);
        
    % Seq 2: Pre-saturation then inversion prior to readout.
    Tis(TisNo).MSatInv = bloch_sim(Pfull, Gfull, RFsatinv, OffResfull, tfull, Tis(TisNo).T1, Tis(TisNo).T2, Ms, LHR, Spoilfull);
    
end

%% Plot
tplot = tfull - tfull(1);
LargeFigWindow(0.5,0.5); 

for SatInv = [false true]
    hold on; col = 'rgbkc';
    clear leg
    for TisNo = 1:3 % Just tissue for now, blood plotted later
        
        if SatInv
            plot(tplot,Tis(TisNo).MSatInv(3,:)   , [col(TisNo) '-'], 'linewidth', 2);
        else
            plot(tplot,Tis(TisNo).MSat(3,:)   , [col(TisNo) '--'], 'linewidth', 2);
        end
        leg{TisNo} = [Tis(TisNo).Name];
    end
    
    % Plot a line where the inversion occurs
    plot([tau tau], [-1 1],'k--');
    
    % Plot a line at zero Mz
    plot([0 max(tplot)],[0 0],'k-')
    
    xlabel 'Time/s'; ylabel 'Mz'; ylim([-1 1]);
    set(gcf,'Position',[521   270   846   356]);
    
end

toLegend(6:8,leg, 'best'); grid on

%% Calculate mean Mz across the readout
LargeFigWindow(0.5,0.2); hold on; col = 'rgbkc'; 
Idx = (tfull >= 0) & (tfull <= 2.0);

BW = 0.2; % Plot bar width

clear leg

for TisNo = 1:3 % Just the tissue components for now (blood plotted later)
    x = 1+(TisNo-3/2-0.5)*BW;
    y = mean(Tis(TisNo).MSat(3,Idx));
    e = std(Tis(TisNo).MSat(3,Idx));
    b = bar(x,y,'facecolor', col(TisNo), 'barwidth',BW);
    hatchfill2(b, 'single', 'HatchAngle', 45, 'HatchDensity', 50, 'HatchColor', 'w', 'HatchLineWidth', 1);

    hold on
    errorbar(x,y,0,e,'color','k')

    x = 2+(TisNo-3/2-0.5)*BW;
    y = mean(Tis(TisNo).MSatInv(3,Idx));
    e = std(Tis(TisNo).MSatInv(3,Idx));
    bar(x,y,'facecolor', col(TisNo), 'barwidth',BW);
    hold on
    errorbar(x,y,0,e,'color','k')
    leg{TisNo} = Tis(TisNo).Name;
end

xlim([0.5 2.5]); ylim([-1 1])
ylabel('Mean Mz')

set(gca,'xtick',1:2)
set(gca,'xticklabel',{'Sat','Sat+Inv'})
xlim([0.5 2.5])
toLegend(1:5:11,leg,'EastOutside'); grid on

%% Plot the blood signal
tplot = tfull - tfull(1);
LargeFigWindow(0.5,0.5); 
dp = 11; % For some reason Matlab rendering dashed lines fails when there are too many data points, so only plot every dp points

for SatInv = [false true] % Plot regular CAPRIA and CAPRIA+S
    hold on; cols = [1 0 0; 0 1 0; 0 0 1; tofadecolor([1 0 0],0.5); 0.5 0 0]; % Use faded and dark red for control/label blood
    clear leg
    for TisNo = 4:5 % Just the blood this time
        
        if SatInv
            plot(tplot(1:dp:end),Tis(TisNo).MSatInv(3,1:dp:end),'-', 'color',cols(TisNo,:), 'linewidth', 2);
        else
            plot(tplot(1:dp:end),Tis(TisNo).MSat(3,1:dp:end)   , '--', 'color',cols(TisNo,:), 'linewidth', 2);
        end
        leg{TisNo} = [Tis(TisNo).Name];
    end
    
    % Plot a line where the inversion occurs
    plot([tau tau], [-1 1],'k--');
    
    % Plot a line at zero Mz
    plot([0 max(tplot)],[0 0],'k-')
    
end

xlabel 'Time/s'; ylabel 'Mz'; ylim([-1 1]);
set(gcf,'Position',[521   270   846   356]);
    
toLegend(5:6,{'Blood (control)','Blood (label)'}, 'best'); grid on