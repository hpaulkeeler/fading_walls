% Reproduce Figures 8-12 from WallsFading manuscript
% Random location model with empirical PDF estimation

clearvars; clc; close all;

% Set random seed for reproducibility
rng(42);

fprintf('Generating Figures 9-13...\n');

%% Parameters
d = 1.0;          % a = b = 0.5, so d = 1.0
beta = 4.0;       % Path loss exponent
kappa = 0.5;      % Absorption coefficient
N_samples = 50000; % Number of Monte Carlo samples
N_plot = 500;     % Number of points for plotting

%% Helper function: Lerch transcendent
function y = lerch_transcendent(zeta, s, alpha, max_terms)
    if nargin < 4
        max_terms = 300;
    end
    
    % Start from n=0 (standard definition)
    y = 0;
    for n = 0:(max_terms-1)
        term = (zeta^n) / ((n + alpha)^s);
        y = y + term;
        if abs(term) < 1e-12
            break;
        end
    end
end

%% Helper function: NLOS signal using Lerch function
function S = nlos_signal_lerch(r, d, kappa, k, beta, max_terms)
    if nargin < 6
        max_terms = 200;
    end
    
    % First term (right images)
    zeta1 = -sqrt(kappa) * exp(1i * k * d);
    alpha1 = -r / d;
    phi1 = lerch_transcendent(zeta1, beta/2, alpha1, max_terms);
    term1 = (exp(-1i * k * r) / (d^(beta/2))) * (phi1 - (-d/r)^(beta/2));
    
    % Second term (left images)
    zeta2 = -sqrt(kappa) * exp(1i * k * d);
    alpha2 = r / d;
    phi2 = lerch_transcendent(zeta2, beta/2, alpha2, max_terms);
    term2 = (exp(1i * k * r) / (d^(beta/2))) * (phi2 - (d/r)^(beta/2));
    
    S = term1 + term2;
end

%% ========================================================================
%% FIGURE 9: Random location in X direction, k=100
%% ========================================================================

fprintf('\nGenerating Figure 9 (k=100, vary x)...\n');

k = 100;

% Top plot: Power function
x_range = linspace(0.15, 0.35, N_plot);
P_power = zeros(size(x_range));

for i = 1:length(x_range)
    r = abs(x_range(i));
    if r < 0.001
        r = 0.001;
    end
    S = nlos_signal_lerch(r, d, kappa, k, beta);
    P_power(i) = abs(S)^2;
end

% Find turning points
turning_idx = islocalmin(P_power) | islocalmax(P_power);
turning_powers = P_power(turning_idx);

% Bottom plot: Empirical PDF via Monte Carlo
x_min = 0.15; x_max = 0.35;
X_samples = x_min + (x_max - x_min) * rand(N_samples, 1);

fprintf('  Sampling %d random locations...\n', N_samples);
P_samples = zeros(N_samples, 1);
for i = 1:N_samples
    r = abs(X_samples(i));
    if r < 0.001
        r = 0.001;
    end
    S = nlos_signal_lerch(r, d, kappa, k, beta);
    P_samples(i) = abs(S)^2;
    
    if mod(i, 10000) == 0
        fprintf('    Processed %d/%d\n', i, N_samples);
    end
end

% Create figure
figure('Position', [100 100 800 900]);

% Top subplot
subplot(2,1,1);
plot(x_range, P_power, 'b-', 'LineWidth', 2); hold on;
plot(x_range(turning_idx), P_power(turning_idx), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
xlabel('x', 'FontSize', 12);
ylabel('Power value P(x,y)', 'FontSize', 12);
title('Vary x, keep y fixed.', 'FontSize', 12);
grid on;
xlim([x_min x_max]);

% Bottom subplot
subplot(2,1,2);
histogram(P_samples, 100, 'Normalization', 'pdf', 'FaceColor', [0.4 0.6 1], 'EdgeColor', 'none');
hold on;
for tp = turning_powers'
    plot([tp tp], ylim, 'r--', 'LineWidth', 1.5);
end
xlabel('Power value s', 'FontSize', 12);
ylabel('Empirical PDF of P', 'FontSize', 12);
grid on;

saveas(gcf, 'Figure9_matlab.png');
fprintf('  Saved Figure9_matlab.png\n');

%% ========================================================================
%% FIGURE 10: Random location in Y direction, k=1000
%% ========================================================================

fprintf('\nGenerating Figure 10 (k=1000, vary y)...\n');

k = 1000;

% Top plot: Power function
y_range = linspace(0.15, 0.35, N_plot);
P_power = zeros(size(y_range));

for i = 1:length(y_range)
    r = abs(y_range(i));
    if r < 0.001
        r = 0.001;
    end
    S = nlos_signal_lerch(r, d, kappa, k, beta, 300);
    P_power(i) = abs(S)^2;
end

% Find turning points
turning_idx = islocalmin(P_power) | islocalmax(P_power);
turning_powers = P_power(turning_idx);

% Bottom plot: Empirical PDF
y_min = 0.15; y_max = 0.35;
Y_samples = y_min + (y_max - y_min) * rand(N_samples, 1);

fprintf('  Sampling %d random locations...\n', N_samples);
P_samples = zeros(N_samples, 1);
for i = 1:N_samples
    r = abs(Y_samples(i));
    if r < 0.001
        r = 0.001;
    end
    S = nlos_signal_lerch(r, d, kappa, k, beta, 300);
    P_samples(i) = abs(S)^2;
    
    if mod(i, 10000) == 0
        fprintf('    Processed %d/%d\n', i, N_samples);
    end
end

% Create figure
figure('Position', [100 100 800 900]);

% Top subplot
subplot(2,1,1);
plot(y_range, P_power, 'b-', 'LineWidth', 2); hold on;
plot(y_range(turning_idx), P_power(turning_idx), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
xlabel('y', 'FontSize', 12);
ylabel('Power value P(x,y)', 'FontSize', 12);
title('Vary y, keep x fixed.', 'FontSize', 12);
grid on;
xlim([y_min y_max]);

% Bottom subplot
subplot(2,1,2);
histogram(P_samples, 100, 'Normalization', 'pdf', 'FaceColor', [0.4 0.6 1], 'EdgeColor', 'none');
hold on;
for tp = turning_powers'
    plot([tp tp], ylim, 'r--', 'LineWidth', 1.5);
end
xlabel('Power value s', 'FontSize', 12);
ylabel('Empirical PDF of P', 'FontSize', 12);
grid on;

saveas(gcf, 'Figure10_matlab.png');
fprintf('  Saved Figure10_matlab.png\n');

%% ========================================================================
%% FIGURE 11: Random location in Y, wider range, k=100
%% ========================================================================

fprintf('\nGenerating Figure 11 (k=100, wider y range)...\n');

k = 100;

% Top plot
y_range = linspace(-0.5, 0.5, N_plot);
P_power = zeros(size(y_range));

for i = 1:length(y_range)
    r = abs(y_range(i));
    if r < 0.001
        r = 0.001;
    end
    S = nlos_signal_lerch(r, d, kappa, k, beta);
    P_power(i) = abs(S)^2;
end

% Find turning points
turning_idx = islocalmin(P_power) | islocalmax(P_power);
turning_powers = P_power(turning_idx);

% Bottom plot
y_min = -0.5; y_max = 0.5;
Y_samples = y_min + (y_max - y_min) * rand(N_samples, 1);

fprintf('  Sampling %d random locations...\n', N_samples);
P_samples = zeros(N_samples, 1);
for i = 1:N_samples
    r = abs(Y_samples(i));
    if r < 0.001
        r = 0.001;
    end
    S = nlos_signal_lerch(r, d, kappa, k, beta);
    P_samples(i) = abs(S)^2;
    
    if mod(i, 10000) == 0
        fprintf('    Processed %d/%d\n', i, N_samples);
    end
end

% Create figure
figure('Position', [100 100 800 900]);

subplot(2,1,1);
plot(y_range, P_power, 'b-', 'LineWidth', 2); hold on;
plot(y_range(turning_idx), P_power(turning_idx), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
xlabel('y', 'FontSize', 12);
ylabel('Power value P(x,y)', 'FontSize', 12);
title('Vary y, keep x fixed.', 'FontSize', 12);
grid on;
xlim([y_min y_max]);

subplot(2,1,2);
histogram(P_samples, 100, 'Normalization', 'pdf', 'FaceColor', [0.4 0.6 1], 'EdgeColor', 'none');
hold on;
for tp = turning_powers'
    plot([tp tp], ylim, 'r--', 'LineWidth', 1.5);
end
xlabel('Power value s', 'FontSize', 12);
ylabel('Empirical PDF of P', 'FontSize', 12);
grid on;

saveas(gcf, 'Figure11_matlab.png');
fprintf('  Saved Figure11_matlab.png\n');

%% ========================================================================
%% FIGURE 12: Random location in Y, smaller interval, k=100
%% ========================================================================

fprintf('\nGenerating Figure 12 (k=100, narrower y range)...\n');

k = 100;

% Top plot
y_range = linspace(0, 0.6, N_plot);
P_power = zeros(size(y_range));

for i = 1:length(y_range)
    r = abs(y_range(i));
    if r < 0.001
        r = 0.001;
    end
    S = nlos_signal_lerch(r, d, kappa, k, beta);
    P_power(i) = abs(S)^2;
end

% Find turning points
turning_idx = islocalmin(P_power) | islocalmax(P_power);
turning_powers = P_power(turning_idx);

% Bottom plot
y_min = 0; y_max = 0.6;
Y_samples = y_min + (y_max - y_min) * rand(N_samples, 1);

fprintf('  Sampling %d random locations...\n', N_samples);
P_samples = zeros(N_samples, 1);
for i = 1:N_samples
    r = abs(Y_samples(i));
    if r < 0.001
        r = 0.001;
    end
    S = nlos_signal_lerch(r, d, kappa, k, beta);
    P_samples(i) = abs(S)^2;
    
    if mod(i, 10000) == 0
        fprintf('    Processed %d/%d\n', i, N_samples);
    end
end

% Create figure
figure('Position', [100 100 800 900]);

subplot(2,1,1);
plot(y_range, P_power, 'b-', 'LineWidth', 2); hold on;
plot(y_range(turning_idx), P_power(turning_idx), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
xlabel('y', 'FontSize', 12);
ylabel('Power value P(x,y)', 'FontSize', 12);
title('Vary y, keep x fixed.', 'FontSize', 12);
grid on;
xlim([y_min y_max]);

subplot(2,1,2);
histogram(P_samples, 100, 'Normalization', 'pdf', 'FaceColor', [0.4 0.6 1], 'EdgeColor', 'none');
hold on;
for tp = turning_powers'
    plot([tp tp], ylim, 'r--', 'LineWidth', 1.5);
end
xlabel('Power value s', 'FontSize', 12);
ylabel('Empirical PDF of P', 'FontSize', 12);
grid on;

saveas(gcf, 'Figure12_matlab.png');
fprintf('  Saved Figure12_matlab.png\n');

%% ========================================================================
%% FIGURE 13: Random location in both X and Y, k=10
%% ========================================================================

fprintf('\nGenerating Figure 13 (k=10, vary x and y)...\n');

k = 10;

% Top plot: 3D surface
N_mesh = 40;
x_vals = linspace(0.05, 0.45, N_mesh);
y_vals = linspace(-0.2, 0.2, N_mesh);
[X_grid, Y_grid] = meshgrid(x_vals, y_vals);
P_surface = zeros(size(X_grid));

fprintf('  Computing 3D surface...\n');
for i = 1:size(X_grid, 1)
    for j = 1:size(X_grid, 2)
        r = sqrt(X_grid(i,j)^2 + Y_grid(i,j)^2);
        if r < 0.001
            r = 0.001;
        end
        S = nlos_signal_lerch(r, d, kappa, k, beta);
        P_surface(i,j) = abs(S)^2;
    end
end

% Bottom plot: Random sampling in 2D
x_min = 0.05; x_max = 0.45;
y_min = -0.2; y_max = 0.2;
X_samples = x_min + (x_max - x_min) * rand(N_samples, 1);
Y_samples = y_min + (y_max - y_min) * rand(N_samples, 1);

fprintf('  Sampling %d random 2D locations...\n', N_samples);
P_samples = zeros(N_samples, 1);
for i = 1:N_samples
    r = sqrt(X_samples(i)^2 + Y_samples(i)^2);
    if r < 0.001
        r = 0.001;
    end
    S = nlos_signal_lerch(r, d, kappa, k, beta);
    P_samples(i) = abs(S)^2;
    
    if mod(i, 10000) == 0
        fprintf('    Processed %d/%d\n', i, N_samples);
    end
end

% Create figure
figure('Position', [100 100 1000 900]);

% Top subplot: 3D surface
subplot(2,1,1);
surf(X_grid, Y_grid, P_surface);
shading interp;
xlabel('x', 'FontSize', 11);
ylabel('y', 'FontSize', 11);
zlabel('P(x,y)', 'FontSize', 11);
title('Vary x and y', 'FontSize', 12);
colorbar;

% Bottom subplot: Empirical PDF
subplot(2,1,2);
histogram(P_samples, 100, 'Normalization', 'pdf', 'FaceColor', [0.4 0.6 1], 'EdgeColor', 'none');
xlabel('s', 'FontSize', 12);
ylabel('Empirical PDF of P', 'FontSize', 12);
title('Random location', 'FontSize', 12);
grid on;

saveas(gcf, 'Figure13_matlab.png');
fprintf('  Saved Figure13_matlab.png\n');

fprintf('\nAll figures generated successfully!\n');
fprintf('\nGenerated files:\n');
fprintf('  - Figure9_matlab.png   (X random, k=100)\n');
fprintf('  - Figure10_matlab.png  (Y random, k=1000)\n');
fprintf('  - Figure11_matlab.png  (Y random, wide range, k=100)\n');
fprintf('  - Figure12_matlab.png  (Y random, narrow range, k=100)\n');
fprintf('  - Figure13_matlab.png  (X,Y both random, k=10)\n');
