% RandomRays.m
% 
% Calculates the sums of reflected electromagnetic rays between two parallel
% walls using the method of images. Generates figures showing both deterministic
% power curves and empirical probability density functions (PDFs) from random
% transmitter locations.
%
% This script implements the random location model from the manuscript 
% "Reflected wireless signals under random spatial sampling" to demonstrate how turning points
% in the deterministic signal power create peaks/singularities in the empirical
% PDF (Proposition VI.1).
%
% Usage:
%   Set choiceVaryXY = 1 for varying x (manuscript Figure 8)
%   Set choiceVaryXY = 2 for varying y (manuscript Figures 9-11)
%   Set choiceVaryXY = 3 for varying x and y (manuscript Figure 12)

clearvars; clc; close all;

% Set default font sizes for publication quality
set(groot, 'DefaultAxesFontSize', 14);
set(groot, 'DefaultTextFontSize', 14);
set(groot, 'DefaultAxesLabelFontSizeMultiplier', 1.15);
set(groot, 'DefaultAxesTitleFontSizeMultiplier', 1.15);
set(groot, 'DefaultLegendFontSize', 12);

%% ========================================================================
%% CONTROL PARAMETERS
%% ========================================================================

% Model selection
choiceVaryXY = 3;       % 1: vary x; 2: vary y; 3: vary x and y (slower)
boolePlotOptima = 1;    % 1 to plot turning points (local extrema); 0 to skip
booleLineSight = 0;     % 1 to include line-of-sight term; 0 for NLOS only
booleDensBins = 1;      % Method for PDF estimation (see funPDF.m)

% Simulation parameters
numbRand = 10^5;        % Number of random samples for empirical PDF
numbBins = 200;         % Number of histogram bins (should be < numbRand)

% Random seed for reproducibility
rng(42);                 % Use rng(42) for manuscript figures

%% ========================================================================
%% PHYSICAL MODEL PARAMETERS
%% ========================================================================

% Wave parameters (one for each choice of choiceVaryXY)
waveNumbValues = [100, 1000, 50];  % Wave number k = 2π/λ
betaPath = 4;                       % Attenuation exponent β (manuscript uses β=4)
kappaAbsorp = 0.5;                  % Absorption coefficient κ ∈ [0,1]

% Wall geometry
a = 0.5;  % Distance from origin to right wall
b = 0.5;  % Distance from origin to left wall

% Transmitter positions (one set for each choiceVaryXY)
xTX = [0.25, 0, 0.25];          % x-coordinate of transmitter center
yTX = [0, 0.25, 0.25];          % y-coordinate of transmitter center
xDeltaTX = [0.2, 0, 0.2];       % x sampling width (±xDelta/2 around xTX)
yDeltaTX = [0, 0.2, 0.2];       % y sampling width (±yDelta/2 around yTX)

% Summation parameters
numbTerm = 100;         % Maximum number of reflection terms to compute
indexFirst = 1;         % Start summation at first reflection (n=1)
indexLast = numbTerm;   % End summation at numbTerm

% Plotting parameters
expPower = 2;           % Exponent for power calculation (P = |S|^expPower)

%% ========================================================================
%% DERIVED PARAMETERS
%% ========================================================================

% Select wave number based on chosen scenario
waveNumb = waveNumbValues(choiceVaryXY);
waveLength = 2*pi/waveNumb;         % Wavelength λ
waveFrequency = 3e8/waveLength;     % Frequency (assuming speed of light)

% Number of points for deterministic plotting
numbPlot = round(1000*waveNumb);

%% ========================================================================
%% GENERATE TRANSMITTER LOCATIONS
%% ========================================================================

% Deterministic locations for power curve plotting
[xPlot, yPlot] = funPosTXPlot(xTX, yTX, xDeltaTX, yDeltaTX, ...
                               choiceVaryXY, numbPlot);

% Random locations for empirical PDF estimation
[xRand, yRand] = funPosTXRand(xTX, yTX, xDeltaTX, yDeltaTX, ...
                               choiceVaryXY, numbRand);

%set distances from origin
if choiceVaryXY==1
    rPlot=xPlot;
    rRand=xRand;
end
if choiceVaryXY==2
    rPlot=yPlot;
    rRand=yRand;
end
if choiceVaryXY==3
    rPlot=hypot(xPlot, yPlot);
    rRand=hypot(xRand, yRand);
end

%% ========================================================================
%% LINE-OF-SIGHT (LOS) SIGNAL
%% ========================================================================

% LOS signal function: S_LOS(r) = r^(-β/2) * exp(j*k*r)
funLOS = @(r0) r0.^(-betaPath/2) .* exp(1i*waveNumb*r0);

% Transmitter location for this scenario
x = xTX(choiceVaryXY);
y = yTX(choiceVaryXY);
r = hypot(x, y);

% Compute LOS term (or set to zero if NLOS-only)
losDet = funLOS(r);

if booleLineSight
    % Include LOS term in both models
    losLoc = losDet;
    losPhase = losDet;
    labelPDF = 'Empirical PDF of $\hat{P}$';      % Full signal (LOS+NLOS)
    labelPower = 'Power value $\hat{P}(x,y)$';
else
    % NLOS only (manuscript uses this setting)
    losDet = 0;
    losLoc = 0;
    losPhase = 0;
    labelPDF = 'Empirical PDF of $P$';            % NLOS signal only
    labelPower = 'Power value $P(x,y)$';
end

%% ========================================================================
%% LABELS FOR PLOTTING
%% ========================================================================

labelTitleCell = {'Vary x, keep y fixed.', ...
                  'Vary y, keep x fixed.', ...
                  'Vary x and y'};
labelXCell = {'$x$', '$y$', '$x$'};
labelYCell = {labelPower, labelPower, '$y$'};
labelZCell = {'', '', labelPower};

% Select labels for current scenario
labelTitle = labelTitleCell{choiceVaryXY};
labelX = labelXCell{choiceVaryXY};
labelY = labelYCell{choiceVaryXY};
labelZ = labelZCell{choiceVaryXY};

%% ========================================================================
%% COMPUTE DETERMINISTIC SIGNAL (for power curve)
%% ========================================================================

% Calculate NLOS signal using method of images
% See funSumNLOS.m for details on the summation
[nlosDet, ~, ~, ~] = funSumNLOS(a, b, xPlot, yPlot, betaPath, ...
                                 waveNumb, kappaAbsorp, indexFirst, indexLast);

% Add LOS term (if included)
sigDet = nlosDet + losDet;

% Compute signal power P(x,y) = |S(x,y)|^2
sigDetPower = abs(sigDet(:)).^expPower;

%% ========================================================================
%% FIGURE 1: Deterministic Power Curve
%% ========================================================================

figure('Name', 'Deterministic Power');

if choiceVaryXY < 3
    % 1D case: Plot power vs position
    set(gcf, 'DefaultLineLineWidth', 2);
    plot(rPlot, sigDetPower);
    grid on;
else
    % 2D case: Plot power surface
    sigDetPower = reshape(sigDetPower, size(rPlot));
    surf(xPlot, yPlot, sigDetPower);
    shading interp;
    zlabel(labelZ, 'Interpreter', 'Latex');
end

xlabel(labelX, 'Interpreter', 'Latex');
ylabel(labelY, 'Interpreter', 'Latex');
title(labelTitle);
axis tight;

%% ========================================================================
%% COMPUTE RANDOM SIGNALS (for empirical PDF)
%% ========================================================================

% Calculate NLOS signals for random transmitter locations
% Returns both random location and random phase models
[nlosLoc, nlosPhase, ~, ~] = funSumNLOS(a, b, xRand, yRand, betaPath, ...
                                         waveNumb, kappaAbsorp, ...
                                         indexFirst, indexLast);

% Add LOS term (if included)
sigLoc = nlosLoc + losLoc;      % Random location model
sigPhase = nlosPhase + losPhase; % Random phase model

% Compute power and separate real/imaginary parts
[sigLocPower, sigLocReal, sigLocImag] = funSignalPower(sigLoc, expPower);
[sigPhasePower, sigPhaseReal, sigPhaseImag] = funSignalPower(sigPhase, expPower);

% Calculate mean power values
mean_sigLocPower = mean(sigLocPower);
mean_sigPhasePower = mean(sigPhasePower);

%% ========================================================================
%% FIGURE 2: Comparison of Random Location and Random Phase Models
%% ========================================================================

figure('Name', 'Random Location vs Random Phase');

axis_xmin = 1.1 * min([sigPhasePower(:); sigLocPower(:)]);
axis_xmax = 0.9 * max([sigPhasePower(:); sigLocPower(:)]);

% Subplot 1: Random location model
subplot(2,1,1);
histLocPower = histogram(sigLocPower, numbBins, 'Normalization', 'pdf');
hold on;
plot(mean_sigLocPower, 0, 'r*', 'MarkerSize', 10);  % Mark mean
xlabel('Power value $s$', 'Interpreter', 'Latex');
ylabel(labelPDF, 'Interpreter', 'Latex');
title('Random location');
grid on;
xlim([axis_xmin axis_xmax]);
axix_ymax = 0.9 * max(histLocPower.Values);
ylim([0 axix_ymax]);

% Subplot 2: Random phase model
subplot(2,1,2);
set(gca, 'FontSize', 14);
histPhasePower = histogram(sigPhasePower, numbBins, 'Normalization', 'pdf');
hold on;
plot(mean_sigPhasePower, 0, 'r*', 'MarkerSize', 10);  % Mark mean
xlabel('Power value $s$', 'Interpreter', 'Latex');
ylabel(labelPDF, 'Interpreter', 'Latex');
title('Random phase');
grid on;
xlim([axis_xmin axis_xmax]);
ylim([0 axix_ymax]);

%% ========================================================================
%% FIGURE 3: Power Curve with Empirical PDF (Main Result Figure)
%% ========================================================================
% This figure demonstrates Proposition VI.1: turning points in the 
% deterministic power function create peaks/singularities in the empirical PDF

figure('Name', 'Power Curve with Turning Points and PDF');

% Subplot 1: Deterministic power curve with turning points marked
subplot(2,1,1);
set(gcf, 'DefaultLineLineWidth', 2);

if choiceVaryXY < 3
    % 1D case: Plot power curve
    plot(rPlot, sigDetPower);
    hold on;
    grid on;
else
    % 2D case: Plot power surface
    surf(xPlot, yPlot, sigDetPower);
    shading interp;
    hold on;
    zlabel(labelPower, 'Interpreter', 'Latex');
end

xlabel(labelX, 'Interpreter', 'Latex');
ylabel(labelY, 'Interpreter', 'Latex');
title(labelTitle);
axis tight;

% Find and plot turning points (local minima and maxima)
if boolePlotOptima
    if choiceVaryXY < 3
        % 1D case: Use islocalmin/islocalmax
        % Note: If Signal Processing Toolbox is available, you can also use:
        %   [pks, locs] = findpeaks(sigDetPower);
        %   [pks_inv, locs_inv] = findpeaks(-sigDetPower);
        
        booleMin = islocalmin(sigDetPower);
        booleMax = islocalmax(sigDetPower);
        booleOptima = booleMin | booleMax;
        optimaLocPower = sigDetPower(booleOptima);
        
        % Plot turning points as red circles
        scatter(rPlot(booleOptima), optimaLocPower, 'r', 'filled');
        
    else
        % 2D case: Use imregionalmin/imregionalmax
        % Note: Requires Image Processing Toolbox
        booleMax = imregionalmax(sigDetPower);
        booleMin = imregionalmin(sigDetPower);
        booleOptima = booleMin | booleMax;
        optimaLocPower = sigDetPower(booleOptima);
        
        % Plot turning points as red circles
        scatter3(xPlot(booleOptima), yPlot(booleOptima), ...
                 optimaLocPower, 'r', 'filled');
    end
end

% Subplot 2: Empirical PDF with turning point values marked
subplot(2,1,2);
histogram(sigLocPower, numbBins, 'Normalization', 'pdf');
hold on;

% Plot vertical lines at turning point power values
if boolePlotOptima
    for tp = optimaLocPower'
        plot([tp tp], ylim, 'r--', 'LineWidth', 1.5);
    end
end

xlabel('Power value $s$', 'Interpreter', 'Latex');
ylabel(labelPDF, 'Interpreter', 'Latex');
title('Empirical PDF from random location sampling');
grid on;

%% ========================================================================
%% DISPLAY SUMMARY STATISTICS
%% ========================================================================

fprintf('\n========== SUMMARY STATISTICS ==========\n');
fprintf('Model parameters:\n');
fprintf('  Wave number k = %g\n', waveNumb);
fprintf('  Path loss exponent β = %g\n', betaPath);
fprintf('  Absorption coefficient κ = %g\n', kappaAbsorp);
fprintf('  Wall distances: a = %g, b = %g\n', a, b);
fprintf('\nDeterministic power curve:\n');
fprintf('  Min power = %.6f\n', min(sigDetPower));
fprintf('  Max power = %.6f\n', max(sigDetPower));
fprintf('  Mean power = %.6f\n', mean(sigDetPower));

if boolePlotOptima && exist('optimaLocPower', 'var')
    fprintf('\nTurning points (N = %d):\n', length(optimaLocPower));
    for i = 1:length(optimaLocPower)
        fprintf('  %2d: %.6f\n', i, optimaLocPower(i));
    end
end

fprintf('\nRandom location model:\n');
fprintf('  Number of samples = %d\n', numbRand);
fprintf('  Mean power = %.6f\n', mean_sigLocPower);
fprintf('  Std power = %.6f\n', std(sigLocPower));
fprintf('  Min power = %.6f\n', min(sigLocPower));
fprintf('  Max power = %.6f\n', max(sigLocPower));

fprintf('\nRandom phase model:\n');
fprintf('  Mean power = %.6f\n', mean_sigPhasePower);
fprintf('  Std power = %.6f\n', std(sigPhasePower));
fprintf('========================================\n\n');


% Assume you've already created/plot into three figures (figure(1), figure(2), figure(3))

% Gather the open figures in a predictable order: 1,2,3 left-to-right
figs = flipud(findobj('Type','figure'));   % flip so Figure 1 is leftmost
if numel(figs) > 3, figs = figs(1:3); end  % only arrange the first three

% Monitor geometry
mp   = get(groot,'MonitorPositions');      % [x y width height] per monitor
rect = mp(1,:);                            % primary monitor
X0 = rect(1); Y0 = rect(2); W = rect(3); H = rect(4);

% Layout parameters
gap = 12;                                  % horizontal gap between windows (pixels)
y   = Y0 + 0;                              % bottom alignment; add margin if needed (e.g., +40)

% Place each figure using its current size
x = X0;
for i = 1:numel(figs)
    set(figs(i),'Units','pixels');         % ensure positions are in pixels
    pos = get(figs(i),'Position');         % [x y w h] current size
    w   = pos(3); h = pos(4);              % preserve width & height
    set(figs(i),'Position',[x, y, w, h]);  % only move: do not resize
    x = x + w + gap;                       % next figure starts to the right
end

