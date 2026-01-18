# Reflected Wireless Signals - MATLAB Code

MATLAB implementation for the manuscript "Reflected wireless signals under random spatial sampling" by H. Paul Keeler.

## Description

This code reproduces the figures and numerical results from the manuscript, which analyzes how randomly positioning transmitters between parallel reflecting walls creates singularities in the probability density function of received signal power.

## Requirements

### MATLAB Version
- MATLAB R2018b or later (tested on R2024a)

### Required Toolboxes
- **Signal Processing Toolbox** - for `islocalmin`, `islocalmax` (finding turning points)
- **Image Processing Toolbox** - for `imregionalmin`, `imregionalmax` (2D turning points in RandomRays.m)

**Note:** The main figure generation script `walls_fading_figures.m` only requires the Signal Processing Toolbox.

## Files

### Main Scripts
- **`walls_fading_figures.m`** - Generates manuscript Figures 8-12 (main results)
- **`RandomRays.m`** - Interactive demonstration of random location model
- **`CompareLerchVsSummation.m`** - Validates Lerch function against direct summation

### Helper Functions
- **`funTruncLerch.m`** - Computes modified Lerch transcendent function
- **`funSumNLOS.m`** - Direct summation using method of images
- **`funPosTXPlot.m`** - Generates deterministic transmitter locations
- **`funPosTXRand.m`** - Generates random transmitter locations
- **`funSignalPower.m`** - Computes signal power and components
- **`funSignalStats.m`** - Computes signal statistics
- **`funPDF.m`** - Estimates probability density functions

## Quick Start

### Generate Manuscript Figures
```matlab
% Generate all main figures (8-12)
walls_fading_figures

% Output: Figure9_matlab.png through Figure13_matlab.png
% Note: Internal numbering differs from manuscript
%   MATLAB Figure 9  = Manuscript Figure 8
%   MATLAB Figure 10 = Manuscript Figure 9
%   MATLAB Figure 11 = Manuscript Figure 10
%   MATLAB Figure 12 = Manuscript Figure 11
%   MATLAB Figure 13 = Manuscript Figure 12
```

### Validate Lerch Function Implementation
```matlab
% Compare Lerch function vs. direct summation
CompareLerchVsSummation

% Should show machine-precision agreement (< 1e-6)
```

### Interactive Demonstration
```matlab
% Explore random location model
RandomRays

% Set choiceVaryXY to control which direction varies:
%   choiceVaryXY = 1: vary x (manuscript Figure 8)
%   choiceVaryXY = 2: vary y (manuscript Figures 9-11)
%   choiceVaryXY = 3: vary x and y (manuscript Figure 12)
```

## Key Parameters

All scripts use these default parameters (matching the manuscript):

```matlab
d = 1.0;          % Wall separation (a = b = 0.5)
beta = 4.0;       % Path loss exponent
kappa = 0.5;      % Absorption coefficient
k = varies;       % Wave number (10, 100, or 1000 depending on figure)
N_samples = 50000; % Monte Carlo samples for PDF estimation
```

## Usage Examples

### Example 1: Generate a Single Figure
```matlab
% Just Figure 8 (X random, k=100)
clearvars; clc; close all;
rng(42);

d = 1.0; beta = 4.0; kappa = 0.5; k = 100;
x_range = linspace(0.15, 0.35, 500);

for i = 1:length(x_range)
    r = abs(x_range(i));
    S = nlos_signal_lerch(r, d, kappa, k, beta);
    P_power(i) = abs(S)^2;
end

plot(x_range, P_power);
xlabel('x'); ylabel('Power P(x,y)');
title('Vary x, keep y fixed');
```

### Example 2: Modify Parameters
```matlab
% Try different wave numbers
RandomRays  % Edit line 51 to set waveNumbValues = [50, 500, 25]
```

### Example 3: Use Helper Functions Directly
```matlab
% Generate random locations and compute NLOS signal
a = 0.5; b = 0.5;
x = 0.2 + 0.1*randn(1000,1);  % Random x-coordinates
y = zeros(1000,1);             % Fixed y
beta = 4.0; k = 100; kappa = 0.5;

nlos = funSumNLOS(a, b, x, y, beta, k, kappa, 1, 200);
P = abs(nlos).^2;

histogram(P, 50, 'Normalization', 'pdf');
xlabel('Power'); ylabel('Empirical PDF');
```

## Implementation Notes

### Lerch Transcendent Function
The code uses the modified Lerch transcendent Φ₁(ζ,s,α) starting from n=1:

```matlab
Φ₁(ζ,s,α) = Σ(n=1 to ∞) ζⁿ/(n+α)ˢ
```

This excludes the n=0 term (LOS signal), which is handled separately.

### Method of Images
The `funSumNLOS.m` implements the classical method of images with:
- 4-case reflection pattern (right-even, left-even, left-odd, right-odd)
- Early termination when terms < 1e-8
- Absorption factor (-√κ)ⁿ per reflection

### Numerical Convergence
- Most cases: 100-200 reflection terms sufficient
- High wave numbers (k > 1000): Use 300-400 terms
- Validation script achieves machine precision (< 1e-12)

## Troubleshooting

### "Undefined function 'islocalmin'"
Install Signal Processing Toolbox:
```matlab
% Check if installed
ver

% If missing, install via Add-Ons in MATLAB GUI
```

### "Undefined function 'imregionalmin'"
Required only for RandomRays.m with 2D mode (choiceVaryXY = 3). Install Image Processing Toolbox or use 1D modes.

### Figures look different from manuscript
Ensure you're using the corrected code with:
- β = 4.0 (not 2.0)
- Correct y-intervals for Figures 10-11
- Random seed: `rng(42)`

## Performance

Approximate runtimes (Intel i7, 16GB RAM):

| Script | Runtime | Memory |
|--------|---------|--------|
| walls_fading_figures.m | ~5-10 min | ~500 MB |
| CompareLerchVsSummation.m | ~5 sec | ~50 MB |
| RandomRays.m (1D) | ~30 sec | ~200 MB |
| RandomRays.m (2D) | ~2 min | ~500 MB |

## Citation

If you use this code, please cite:

```bibtex
@article{keeler2026walls,
  title={Reflected wireless signals under random spatial sampling},
  author={Keeler, H. Paul},
  journal={Preprint},
  year={2026}
}
```

## License

MIT License - see LICENSE file for details

## Acknowledgments

This work was developed with assistance from Claude (Anthropic) for code documentation and organization.
