% CompareLerchVsSummation.m
%
% Demonstrates mathematical equivalence between two methods for computing
% NLOS (non-line-of-sight) electromagnetic signals between parallel walls:
%   1. Direct summation using method of images
%   2. Closed-form expression using modified Lerch transcendent function
%
% The modified Lerch function Φ₁(ζ,s,α) is defined as:
%
%   Φ₁(ζ,s,α) = Σ_{n=1}^∞ ζ^n / (n + α)^s
%
% This differs from the standard Lerch Φ(ζ,s,α) by starting at n=1 instead
% of n=0. The n=0 term represents the direct line-of-sight (LOS) signal,
% which is handled separately, making Φ₁ appropriate for NLOS calculations.
%
% REFERENCE:
%   Manuscript "Reflected wireless signals under random spatial sampling"
%   Proposition V.1 (Equation 29): Closed-form NLOS signal expression

clearvars; clc; close all;

% Set default font sizes for publication quality
set(groot, 'DefaultAxesFontSize', 14);
set(groot, 'DefaultTextFontSize', 14);
set(groot, 'DefaultAxesLabelFontSizeMultiplier', 1.15);
set(groot, 'DefaultAxesTitleFontSizeMultiplier', 1.15);
set(groot, 'DefaultLegendFontSize', 12);

fprintf('========================================\n');
fprintf('Comparing Lerch Function vs Direct Sum\n');
fprintf('========================================\n\n');

%% ========================================================================
%% PARAMETERS (Using manuscript Figure 9 settings)
%% ========================================================================

% Physical parameters
a = 0.5;                % Distance to right wall
b = 0.5;                % Distance to left wall
d = a + b;              % Total wall separation
betaPath = 4;           % Path loss exponent β
kappaAbsorp = 0.5;      % Absorption coefficient κ
waveNumb = 100;         % Wavenumber k

% Summation parameters
indexFirst = 1;         % Start from first reflection (n=1)
indexLast = 200;        % Use 200 terms for convergence

% Test points
N_test = 10;            % Number of test points
x_test = linspace(0.15, 0.35, N_test);  % x coordinates
y_test = zeros(1, N_test);              % y = 0 (fixed)

%% ========================================================================
%% METHOD 1: Direct Summation (Method of Images)
%% ========================================================================

fprintf('METHOD 1: Direct summation using method of images\n');
fprintf('--------------------------------------------------\n');

tic;
[nlosDirect, ~, ~, ~] = funSumNLOS(a, b, x_test, y_test, betaPath, ...
                                    waveNumb, kappaAbsorp, indexFirst, indexLast);
time_direct = toc;

fprintf('Computed %d test points in %.4f seconds\n', N_test, time_direct);
fprintf('Average per point: %.6f seconds\n\n', time_direct/N_test);

%% ========================================================================
%% METHOD 2: Closed-Form Using Modified Lerch Function Φ₁
%% ========================================================================

fprintf('METHOD 2: Closed-form using modified Lerch function Φ₁\n');
fprintf('------------------------------------------------------\n');

% Preallocate
nlosLerch = zeros(N_test, 1);

tic;
for i = 1:N_test
    % Distance from origin
    r = hypot(x_test(i), y_test(i));
    if r < 0.001
        r = 0.001;  % Avoid singularity
    end
    
    % Common parameters for Lerch function
    zeta = -sqrt(kappaAbsorp) * exp(1i * waveNumb * d);
    s = betaPath / 2;
    
    % First term: right images (Equation 29a in manuscript)
    alpha1 = -r / d;
    phi1 = funTruncLerch(zeta, s, alpha1);
    term1 = (exp(-1i * waveNumb * r) / (d^s)) * phi1;
    
    % Second term: left images (Equation 29b in manuscript)
    alpha2 = r / d;
    phi2 = funTruncLerch(zeta, s, alpha2);
    term2 = (exp(1i * waveNumb * r) / (d^s)) * phi2;
    
    % Total NLOS signal
    nlosLerch(i) = term1 + term2;
end
time_lerch = toc;

fprintf('Computed %d test points in %.4f seconds\n', N_test, time_lerch);
fprintf('Average per point: %.6f seconds\n', time_lerch/N_test);
fprintf('Speedup factor: %.2fx faster\n\n', time_direct/time_lerch);

%% ========================================================================
%% COMPARISON: Compute Differences
%% ========================================================================

fprintf('COMPARISON RESULTS\n');
fprintf('--------------------------------------------------\n');

% Compute differences
diff_complex = nlosDirect - nlosLerch;
diff_magnitude = abs(diff_complex);
diff_power = abs(nlosDirect).^2 - abs(nlosLerch).^2;

% Statistical summary
fprintf('Complex signal differences:\n');
fprintf('  Max |S_direct - S_Lerch| = %.2e\n', max(diff_magnitude));
fprintf('  Mean |S_direct - S_Lerch| = %.2e\n', mean(diff_magnitude));
fprintf('  Relative error = %.2e%%\n\n', 100*max(diff_magnitude)/mean(abs(nlosDirect)));

fprintf('Power differences:\n');
fprintf('  Max |P_direct - P_Lerch| = %.2e\n', max(abs(diff_power)));
fprintf('  Mean |P_direct - P_Lerch| = %.2e\n', mean(abs(diff_power)));
fprintf('  Relative error = %.2e%%\n\n', 100*max(abs(diff_power))/mean(abs(nlosDirect).^2));

%% ========================================================================
%% DETAILED COMPARISON TABLE
%% ========================================================================

fprintf('DETAILED POINT-BY-POINT COMPARISON\n');
fprintf('--------------------------------------------------\n');
fprintf('Point |   x   |   Direct Method   |   Lerch Method    | Difference\n');
fprintf('------|-------|-------------------|-------------------|-----------\n');

for i = 1:N_test
    r = hypot(x_test(i), y_test(i));
    mag_direct = abs(nlosDirect(i));
    mag_lerch = abs(nlosLerch(i));
    diff = abs(nlosDirect(i) - nlosLerch(i));
    
    fprintf('%5d | %.3f | %8.6f + %8.6fi | %8.6f + %8.6fi | %.2e\n', ...
            i, x_test(i), real(nlosDirect(i)), imag(nlosDirect(i)), ...
            real(nlosLerch(i)), imag(nlosLerch(i)), diff);
end
fprintf('\n');

%% ========================================================================
%% VISUALIZATION: Signal Magnitude and Power
%% ========================================================================

figure('Position', [100 100 1200 800], 'Name', 'Lerch vs Direct Summation Comparison');

% Compute signal magnitudes and power
mag_direct = abs(nlosDirect);
mag_lerch = abs(nlosLerch);
power_direct = mag_direct.^2;
power_lerch = mag_lerch.^2;

% Subplot 1: Signal magnitude comparison
subplot(2,2,1);
plot(x_test, mag_direct, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Direct Sum'); 
hold on;
plot(x_test, mag_lerch, 'rx--', 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'Lerch Φ₁');
xlabel('x position', 'FontSize', 14);
ylabel('Signal magnitude |S|', 'FontSize', 14);
title('Signal Magnitude: Two Methods', 'FontSize', 16);
legend('Location', 'best', 'FontSize', 12);
grid on;

% Subplot 2: Magnitude difference (zoomed)
subplot(2,2,2);
plot(x_test, diff_magnitude, 'ko-', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('x position', 'FontSize', 14);
ylabel('|S_{direct} - S_{Lerch}|', 'FontSize', 14);
title('Magnitude Difference (max error shown)', 'FontSize', 16);
grid on;
ylim([0, max(diff_magnitude)*1.2]);
text(mean(x_test), max(diff_magnitude)*0.9, ...
     sprintf('Max error: %.2e', max(diff_magnitude)), ...
     'HorizontalAlignment', 'center', 'FontSize', 12, 'BackgroundColor', 'white');

% Subplot 3: Power comparison
subplot(2,2,3);
plot(x_test, power_direct, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Direct Sum');
hold on;
plot(x_test, power_lerch, 'rx--', 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'Lerch Φ₁');
xlabel('x position', 'FontSize', 14);
ylabel('Signal power P = |S|^2', 'FontSize', 14);
title('Signal Power: Two Methods', 'FontSize', 16);
legend('Location', 'best', 'FontSize', 12);
grid on;

% Subplot 4: Power difference
subplot(2,2,4);
plot(x_test, abs(diff_power), 'ko-', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('x position', 'FontSize', 14);
ylabel('|P_{direct} - P_{Lerch}|', 'FontSize', 14);
title('Power Difference', 'FontSize', 16);
grid on;
ylim([0, max(abs(diff_power))*1.2]);
text(mean(x_test), max(abs(diff_power))*0.9, ...
     sprintf('Max error: %.2e', max(abs(diff_power))), ...
     'HorizontalAlignment', 'center', 'FontSize', 12, 'BackgroundColor', 'white');

%% ========================================================================
%% MATHEMATICAL EXPLANATION
%% ========================================================================

fprintf('MATHEMATICAL EQUIVALENCE\n');
fprintf('========================================\n\n');

fprintf('The NLOS signal can be expressed as an infinite sum over reflections:\n\n');
fprintf('  S_NLOS(r) = Σ_{n=1}^∞ (-√κ)^n × ℓ(r̂_n) × e^(jkr̂_n)\n');
fprintf('              + Σ_{n=1}^∞ (-√κ)^n × ℓ(t̂_n) × e^(jkt̂_n)\n\n');

fprintf('where r̂_n and t̂_n are distances to right and left image sources,\n');
fprintf('and ℓ(r) = r^(-β/2) is the path loss function.\n\n');

fprintf('For the symmetric case (a = b = d/2), this simplifies to:\n\n');
fprintf('  S_NLOS(r) = (e^(-jkr) / d^(β/2)) × [Φ₁(ζ, β/2, -r/d)]\n');
fprintf('              + (e^(+jkr) / d^(β/2)) × [Φ₁(ζ, β/2, +r/d)]\n\n');

fprintf('where:\n');
fprintf('  ζ = -√κ × e^(jkd)  (absorption and phase per round-trip)\n');
fprintf('  Φ₁(ζ,s,α) = Σ_{n=1}^∞ ζ^n/(n+α)^s  (modified Lerch function)\n\n');

fprintf('KEY INSIGHT:\n');
fprintf('------------\n');
fprintf('The modified Lerch Φ₁ starts from n=1 (not n=0) because:\n');
fprintf('  • The n=0 term represents the direct LOS path\n');
fprintf('  • NLOS calculations exclude this direct path\n');
fprintf('  • This is mathematically equivalent to using standard Lerch Φ\n');
fprintf('    with explicit subtraction: Φ(ζ,s,α) - α^(-s)\n\n');

fprintf('Convergence: Both methods use %d terms. Increasing to 300-400 terms\n', indexLast);
fprintf('may improve accuracy for large wave numbers (k > 1000).\n\n');

%% ========================================================================
%% CONCLUSION
%% ========================================================================

fprintf('CONCLUSION\n');
fprintf('========================================\n\n');

if max(diff_magnitude) < 1e-6
    fprintf('✓ VALIDATION SUCCESSFUL\n\n');
    fprintf('The two methods produce identical results within numerical precision:\n');
    fprintf('  • Maximum difference: %.2e\n', max(diff_magnitude));
    fprintf('  • Relative error: %.2e%%\n\n', 100*max(diff_magnitude)/mean(mag_direct));
    fprintf('The closed-form Lerch expression (Proposition V.1) is confirmed to be\n');
    fprintf('mathematically equivalent to the direct method-of-images summation.\n');
else
    fprintf('⚠ WARNING: Significant difference detected\n\n');
    fprintf('Consider increasing number of terms (indexLast) for better convergence.\n');
end

fprintf('\n========================================\n');
