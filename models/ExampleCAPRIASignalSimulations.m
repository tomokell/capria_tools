% Example signal simulations for CAPRIA 
%
% Tom Okell, June 2022

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

%% Set sequence parameters
tau = 1.4;          % Labelling duration (s)
CFAAlpha = 6;       % CFA flip angle (degs)
VFAAlpha = [2 9];   % VFA flip angles start and finish (alpha1 and alphaN, degs)
TR = 9e-3;          % TR (ms)
NsegsA = 24;        % Number of segments (number of excitation pulses within temporal resolution) for angio
NphasesA = 9;       % Number of frames for angio
NsegsP = 36;        % Number of segments (number of excitation pulses within temporal resolution) for perfusion
NphasesP = 6;       % Number of frames for perfusion

% Total readout time
T = TR*NsegsA*NphasesA; % NB. Same as TR*NsegsP*NphasesP

% Define the start of imaging here to be right after the labelling finishes
t0 = tau;

% Angiographic physiological parameters
Aparams.T1b = 1.65;     % T1 of blood (s)
Aparams.v = 1;          % Volume fraction occupied by macrovasculature
Aparams.delta_t = 0.75; % Macrovascular transit time (s)

% Perfusion physiological parameters
Pparams.T1b = 1.65;     % T1 of blood (s)
Pparams.Deltat = 1.5;   % Arterial transit time to tissue (s)
Pparams.T1 = 1.3;       % T1 of tissue (s)
Pparams.f = 60;         % CBF (ml/100g/min)

% Start time for simulation
tstart = 0.5; % s

% Angio CFA
[CFA_SigA, t, CFA_Alpha] = CAPRIASignalSimple('Angio','CFA'      ,CFAAlpha,tstart,tau,t0,NsegsA,NphasesA,TR,Aparams,true,true);

% Angio VFA
[VFA_SigA, ~, VFA_Alpha] = CAPRIASignalSimple('Angio','Quadratic',VFAAlpha,tstart,tau,t0,NsegsA,NphasesA,TR,Aparams,true,true);

% Perfusion CFA
CFA_SigP                 = CAPRIASignalSimple('Perf' ,'CFA'      ,CFAAlpha,tstart,tau,t0,NsegsP,NphasesP,TR,Pparams,true,true);

% Perfusion VFA
VFA_SigP                 = CAPRIASignalSimple('Perf' ,'Quadratic',VFAAlpha,tstart,tau,t0,NsegsP,NphasesP,TR,Pparams,true,true);


%% Optimise across a range of VFA values and physiological parameters
% Vary starting and final flip angles:
Alpha1 = 0.5:0.5:20;
AlphaN = 0.5:0.5:30;

% And vary angio and perfusion arrival times (CBF and blood volume just
% scale the signal, so won't feed into the optimisation)
Angio_delta_t = 0.2:0.1:1;
Perf_Delta_t = 0.5:0.1:2;

AngioMeanSig = zeros(length(Alpha1),length(AlphaN),length(Angio_delta_t));
PerfMeanSig  = zeros(length(Alpha1),length(AlphaN),length(Perf_Delta_t ));

% Loop over possible flip angles
for ii = 1:length(Alpha1)
    for jj = 1:length(AlphaN)
        VFAAlpha = [Alpha1(ii) AlphaN(jj)];
        
        % Simulate across all angiographic transit times, outputting only
        % the magnetisation averaged over the imaging time points
        for kk = 1:length(Angio_delta_t)
            Aparams.delta_t = Angio_delta_t(kk);
            [~, ~, ~, ~, ~, tAv, dMAv] = CAPRIASignalSimple('Angio','Quadratic',VFAAlpha,tstart,tau,t0,NsegsA,NphasesA,TR,Aparams,true,true);
            
            % We only want to optimise for signal when it's present in the
            % voxel:
            Idx = (tAv >= Angio_delta_t(kk)) & (tAv < tau + Angio_delta_t(kk)); 
            
            % For very early arrival, we may not hit the criteria above, so
            % make sure at least the maximum signal is included in the
            % index:
            Idx(dMAv==max(dMAv))=true;
            
            % Calculate the mean across these timepoints
            AngioMeanSig(ii,jj,kk) = mean(dMAv(Idx));
        end
        
        % Repeat for perfusion
        for kk = 1:length(Perf_Delta_t)
            Pparams.Deltat = Perf_Delta_t(kk);
            [~, ~, ~, ~, ~, tAv, dMAv] = CAPRIASignalSimple('Perf' ,'Quadratic',VFAAlpha,tstart,tau,t0,NsegsP,NphasesP,TR,Pparams,true,true);            
            
            % Here we take all time points after the peak, where the signal
            % is approximately CBF weighted
            Idx = (tAv >= tau+Perf_Delta_t(kk)); 
            
            % As above, in case of very late arrival, make sure we at least
            % include the final timepoint
            Idx(end) = true;
            
            % Take the average across these timepoints
            PerfMeanSig(ii,jj,kk) = mean(dMAv(Idx) );        
        end
    end
end

%% Take the average across physiological parameters and normalise
% Angiographic metric
AngioMetricForPlot = mean(AngioMeanSig,3); 
AngioMetricForPlot = AngioMetricForPlot/max(AngioMetricForPlot(:));

% Perfusion metric
PerfMetricForPlot = mean(PerfMeanSig,3); 
PerfMetricForPlot = PerfMetricForPlot/max(PerfMetricForPlot(:));

% Combined metric (upweights perfusion, which is lower SNR)
CombMetricForPlot = 0.5*AngioMetricForPlot/max(AngioMetricForPlot(:)) + PerfMetricForPlot/max(PerfMetricForPlot(:));
CombMetricForPlot = CombMetricForPlot/max(CombMetricForPlot(:));


%% Plot the example signal timecourses and flip angle schedules
LargeFigWindow(0.8,0.9); subplot(2,3,1); 
plot(t,CFA_Alpha,t,VFA_Alpha,'linewidth',2); title 'A) Example flip angle schedules'
hold on; ylims = get(gca,'ylim'); plot([t0 t0],ylims,'k--'); xlim([min(t) max(t)])
xlabel 'Time/s'; ylabel 'Flip angle, \alpha/\circ'; xlim([min(t) max(t)])

subplot(2,3,2)
plot(t,CFA_SigA, t, VFA_SigA,'linewidth',2); title 'B) Example angiographic signals'
hold on; ylims = get(gca,'ylim'); plot([t0 t0],ylims,'k--'); xlim([min(t) max(t)])
ylabel 'ASL signal'; xlabel 'Time/s'; 

subplot(2,3,3); 
plot(t,CFA_SigP, t, VFA_SigP,'linewidth',2); title 'C) Example perfusion signals'
hold on; ylims = get(gca,'ylim'); plot([t0 t0],ylims,'k--'); xlim([min(t) max(t)])
ylabel 'ASL signal'; xlabel 'Time/s'; 
leg = {'CFA (6\circ)','VFA (2-9\circ)'}; legend(leg,'location','north')

disp(['Signal amplification at last time point is ' ns(VFA_SigP(end)/CFA_SigP(end))])

% Plot the contours of the optimisation metrics
subplot(2,3,4); contourf(Alpha1,AlphaN,AngioMetricForPlot',15); colorbar; 
hold on; plot(Alpha1,Alpha1,'k--'); 
plot(2,9,'ro','MarkerFaceColor','r'); plot(6,6,'ko','MarkerFaceColor','k'); 
xlabel 'Starting flip angle, \alpha_1/\circ'; ylabel 'Final flip angle, \alpha_N/\circ';
title 'D) Angiography optimisation'; caxis([0 1])

subplot(2,3,5); contourf(Alpha1,AlphaN,PerfMetricForPlot',15); colorbar; 
hold on; plot(Alpha1,Alpha1,'k--'); 
plot(2,9,'ro','MarkerFaceColor','r'); plot(6,6,'ko','MarkerFaceColor','k'); 
xlabel 'Starting flip angle, \alpha_1/\circ'; ylabel 'Final flip angle, \alpha_N/\circ';
title 'E) Perfusion optimisation'; caxis([0 1])

subplot(2,3,6); [C,h] = contourf(Alpha1,AlphaN,CombMetricForPlot',0:0.05:1); colorbar; 
clabel(C,h,0.95,'LabelSpacing',100,'Margin',1);
hold on; plot(Alpha1,Alpha1,'k--'); 
plot(2,9,'ro','MarkerFaceColor','r'); plot(6,6,'ko','MarkerFaceColor','k'); 
xlabel 'Starting flip angle, \alpha_1/\circ'; ylabel 'Final flip angle, \alpha_N/\circ';
title 'F) Combined optimisation'; caxis([0 1])
toLegend([2 3 4],{'CFA constraint','VFA optimum','CFA optimum'},'best');

%% Example signal simulations incorporating dispersion
VFAAlpha = [2 9];   % VFA flip angles start and finish (alpha1 and alphaN, degs)

% Angiographic physiological parameters
Aparams.v = 1;          % Volume fraction occupied by macrovasculature
Aparams.delta_t = 0.75; % Macrovascular transit time (s)
Aparams.s    =   10; % s^-1
Aparams.p    =   0.1; % s

% Perfusion physiological parameters
Pparams.Deltat = 1.5;   % Arterial transit time to tissue (s)
Pparams.f = 60;         % CBF (ml/100g/min)
Pparams.s = Aparams.s;
Pparams.p = Aparams.p;

% Angio VFA with dispersion
VFA_SigAdisp = CAPRIASignalDisp('Angio','Quadratic',VFAAlpha,tstart,tau,t0,NsegsA,NphasesA,TR,Aparams,true,true);

% Perfusion VFA with dispersion
VFA_SigPdisp = CAPRIASignalDisp('Perf' ,'Quadratic',VFAAlpha,tstart,tau,t0,NsegsP,NphasesP,TR,Pparams,true,true);

%% Plot
LargeFigWindow(0.4,0.9); subplot(1,2,1); 
plot(t,VFA_SigA, t, VFA_SigAdisp,'linewidth',2); title 'Example angiographic signals'
hold on; ylims = get(gca,'ylim'); plot([t0 t0],ylims,'k--'); xlim([min(t) max(t)])
ylabel 'ASL signal'; xlabel 'Time/s'; 

subplot(1,2,2)
plot(t,VFA_SigP, t, VFA_SigPdisp,'linewidth',2); title 'Example perfusion signals'
hold on; ylims = get(gca,'ylim'); plot([t0 t0],ylims,'k--'); xlim([min(t) max(t)])
ylabel 'ASL signal'; xlabel 'Time/s'; 
leg = {'No dispersion','With dispersion'}; legend(leg,'location','north')