#!/usr/bin/env python3
"""
RandomRays.py

Calculates the sums of reflected electromagnetic rays between two parallel
walls using the method of images. Generates figures showing both deterministic
power curves and empirical probability density functions (PDFs) from random
transmitter locations.

This script implements the random location model from the manuscript
"Reflected wireless signals under random spatial sampling" to demonstrate how turning points
in the deterministic signal power create peaks/singularities in the empirical
PDF (Proposition VI.1).

Usage:
    Set choiceVaryXY = 0 for varying x (manuscript Figure 8)
    Set choiceVaryXY = 1 for varying y (manuscript Figures 9-11)
    Set choiceVaryXY = 2 for varying x and y (manuscript Figure 12)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks  # Optional: requires scipy
import time

# Set default font sizes for publication quality
plt.rcParams.update({
    'font.size': 14,
    'axes.labelsize': 16,
    'axes.titlesize': 16,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
})

# Import helper functions from separate module
from walls_helpers import (
    funPosTXPlot,
    funPosTXRand,
    funSumNLOS,
    funSignalPower
)

# =========================================================================
# CONTROL PARAMETERS
# =========================================================================

# Model selection
choiceVaryXY = 1       # 0: vary x; 1: vary y; 2: vary x and y (slower)
boolePlotOptima = True  # True to plot turning points; False to skip
booleLineSight = False  # True to include line-of-sight term; False for NLOS only
booleDensBins = True    # Method for PDF estimation (see funPDF)

# Simulation parameters
numbRand = 100000       # Number of random samples for empirical PDF
numbBins = 200          # Number of histogram bins (should be < numbRand)

# Random seed for reproducibility
np.random.seed(42)        # Use seed(3) for manuscript figures

# =========================================================================
# PHYSICAL MODEL PARAMETERS
# =========================================================================

# Wave parameters (one for each choice of choiceVaryXY)
waveNumbValues = np.array([100, 1000, 50])  # Wave number k = 2π/λ
betaPath = 4.0                                # Attenuation exponent β (manuscript uses β=4)
kappaAbsorp = 0.5                             # Absorption coefficient κ ∈ [0,1]

# Wall geometry
a = 0.5  # Distance from origin to right wall
b = 0.5  # Distance from origin to left wall

# Transmitter positions (one set for each choiceVaryXY)
xTX = np.array([0.25, 0, 0.25])          # x-coordinate of transmitter center
yTX = np.array([0, 0.25, 0.25])          # y-coordinate of transmitter center
xDeltaTX = np.array([0.2, 0, 0.2])      # x sampling width (±xDelta/2 around xTX)
yDeltaTX = np.array([0, 0.2, 0.2])      # y sampling width (±yDelta/2 around yTX)

# Summation parameters
numbTerm = 100          # Maximum number of reflection terms to compute
indexFirst = 1          # Start summation at first reflection (n=1)
indexLast = numbTerm   # End summation at numbTerm

# Plotting parameters
expPower = 2            # Exponent for power calculation (P = |S|^expPower)

# =========================================================================
# DERIVED PARAMETERS
# =========================================================================

# Select wave number based on chosen scenario
waveNumb = waveNumbValues[choiceVaryXY]
waveLength = 2*np.pi/waveNumb         # Wavelength λ
waveFrequency = 3e8/waveLength        # Frequency (assuming speed of light)

# Number of points for deterministic plotting
numbPlot = round(200*waveNumb)

# =========================================================================
# GENERATE TRANSMITTER LOCATIONS
# =========================================================================

print("Generating transmitter locations...")

# Deterministic locations for power curve plotting
xPlot, yPlot = funPosTXPlot(xTX, yTX, xDeltaTX, yDeltaTX,
                             choiceVaryXY, numbPlot)

# Random locations for empirical PDF estimation
xRand, yRand = funPosTXRand(xTX, yTX, xDeltaTX, yDeltaTX,
                             choiceVaryXY, numbRand)

# Set distances from origin
if choiceVaryXY == 0:
    rPlot = xPlot
    rRand = xRand
elif choiceVaryXY == 1:
    rPlot = yPlot
    rRand = yRand
elif choiceVaryXY == 2:
    rPlot = np.hypot(xPlot, yPlot)
    rRand = np.hypot(xRand, yRand)

# =========================================================================
# LINE-OF-SIGHT (LOS) SIGNAL
# =========================================================================

# LOS signal function: S_LOS(r) = r^(-β/2) * exp(j*k*r)
funLOS = lambda r0: r0**(-betaPath/2) * np.exp(1j*waveNumb*r0)

# Transmitter location for this scenario
x = xTX[choiceVaryXY]
y = yTX[choiceVaryXY]
r = np.hypot(x, y)

# Compute LOS term (or set to zero if NLOS-only)
losDet = funLOS(r)

if booleLineSight:
    losLoc = losDet
    labelPDF = 'Empirical PDF of $\\hat{P}$'
    labelPower = 'Power value $\\hat{P}(x,y)$'
else:
    losDet = 0
    losLoc = 0
    labelPDF = 'Empirical PDF of $P$'
    labelPower = 'Power value $P(x,y)$'

# =========================================================================
# LABELS FOR PLOTTING
# =========================================================================

labelTitleList = ['Vary x, keep y fixed.',
                  'Vary y, keep x fixed.',
                  'Vary x and y']
labelXList = ['$x$', '$y$', '$x$']
labelYList = [labelPower, labelPower, '$y$']
labelZList = ['', '', labelPower]

labelTitle = labelTitleList[choiceVaryXY]
labelX = labelXList[choiceVaryXY]
labelY = labelYList[choiceVaryXY]
labelZ = labelZList[choiceVaryXY]

# =========================================================================
# COMPUTE DETERMINISTIC SIGNAL (for power curve)
# =========================================================================

print("Computing deterministic signal...")
startTime = time.time()

# Calculate NLOS signal using method of images
nlosDet = funSumNLOS(a, b, xPlot, yPlot, betaPath,
                      waveNumb, kappaAbsorp, indexFirst, indexLast)

# Add LOS term (if included)
sigDet = nlosDet + losDet

# Compute signal power P(x,y) = |S(x,y)|^2
sigDetPower = np.abs(sigDet.flatten())**expPower

print(f"  Computed in {time.time() - startTime:.2f} seconds")

# =========================================================================
# FIGURE 1: Deterministic Power Curve
# =========================================================================

print("Creating Figure 1...")

fig1 = plt.figure(figsize=(10, 6))
fig1.canvas.manager.set_window_title('Deterministic Power')

if choiceVaryXY < 2:
    # 1D case: Plot power vs position
    plt.plot(rPlot, sigDetPower, 'b-', linewidth=2)
    plt.grid(True)
else:
    # 2D case: Plot power surface
    from mpl_toolkits.mplot3d import Axes3D
    ax = fig1.add_subplot(111, projection='3d')
    sigDetPower2D = sigDetPower.reshape(xPlot.shape)
    ax.plot_surface(xPlot, yPlot, sigDetPower2D, cmap='viridis')
    ax.set_zlabel(labelZ)

plt.xlabel(labelX)
plt.ylabel(labelY)
plt.title(labelTitle)
plt.tight_layout()

# =========================================================================
# COMPUTE RANDOM SIGNALS (for empirical PDF)
# =========================================================================

print(f"Computing random signals for {numbRand} samples...")
startTime = time.time()

# Calculate NLOS signals for random transmitter locations
nlosLoc = funSumNLOS(a, b, xRand, yRand, betaPath,
                      waveNumb, kappaAbsorp, indexFirst, indexLast)

# Add LOS term (if included)
sigLoc = nlosLoc + losLoc

# Compute power
sigLocPower, sigLocReal, sigLocImag = funSignalPower(sigLoc, expPower)

# Calculate mean power
meanSigLocPower = np.mean(sigLocPower)

print(f"  Computed in {time.time() - startTime:.2f} seconds")

# =========================================================================
# FIGURE 2: Power Curve with Empirical PDF (Main Result Figure)
# =========================================================================

print("Creating Figure 2...")

fig2 = plt.figure(figsize=(10, 10))
fig2.canvas.manager.set_window_title('Power Curve with Turning Points and PDF')

# Subplot 1: Deterministic power curve with turning points marked
ax1 = plt.subplot(2, 1, 1)

if choiceVaryXY < 2:
    # 1D case: Plot power curve
    plt.plot(rPlot, sigDetPower, 'b-', linewidth=2)
    plt.grid(True)
else:
    # 2D case: Plot power surface
    from mpl_toolkits.mplot3d import Axes3D
    ax1 = fig2.add_subplot(2, 1, 1, projection='3d')
    sigDetPower2D = sigDetPower.reshape(xPlot.shape)
    ax1.plot_surface(xPlot, yPlot, sigDetPower2D, cmap='viridis')

plt.xlabel(labelX, fontsize=16)
plt.ylabel(labelY, fontsize=16)
plt.title(labelTitle, fontsize=16)

# Find and plot turning points (local minima and maxima)
if boolePlotOptima:
    if choiceVaryXY < 2:
        # Method 1: Using scipy.signal.find_peaks (requires scipy)
        # Uncomment if scipy is available:
        # peaksMax, _ = find_peaks(sigDetPower, prominence=0.01)
        # peaksMin, _ = find_peaks(-sigDetPower, prominence=0.01)
        # booleOptima = np.zeros(len(sigDetPower), dtype=bool)
        # booleOptima[peaksMax] = True
        # booleOptima[peaksMin] = True

        # Method 2: Manual detection (always works)
        booleMin = np.zeros(len(sigDetPower), dtype=bool)
        booleMax = np.zeros(len(sigDetPower), dtype=bool)

        for i in range(1, len(sigDetPower)-1):
            if sigDetPower[i] < sigDetPower[i-1] and sigDetPower[i] < sigDetPower[i+1]:
                booleMin[i] = True
            if sigDetPower[i] > sigDetPower[i-1] and sigDetPower[i] > sigDetPower[i+1]:
                booleMax[i] = True

        booleOptima = booleMin | booleMax
        optimaLocPower = sigDetPower[booleOptima]

        # Plot turning points as red circles
        plt.scatter(rPlot[booleOptima], optimaLocPower,
                   c='r', s=50, zorder=5, label='Turning points')

# Subplot 2: Empirical PDF with turning point values marked
ax2 = plt.subplot(2, 1, 2)
histCounts, binEdges, _ = plt.hist(sigLocPower, bins=numbBins,
                                    density=True, alpha=0.7,
                                    color='skyblue', edgecolor='none')

# Plot vertical lines at turning point power values
if boolePlotOptima and choiceVaryXY < 2:
    yLim = plt.ylim()
    for tp in optimaLocPower:
        plt.plot([tp, tp], yLim, 'r--', linewidth=1.5, alpha=0.7)

plt.xlabel('Power value $s$', fontsize=16)
plt.ylabel(labelPDF, fontsize=16)
plt.title('Empirical PDF from random location sampling', fontsize=16)
plt.grid(True, alpha=0.3)

plt.tight_layout()

# =========================================================================
# DISPLAY SUMMARY STATISTICS
# =========================================================================

print('\n========== SUMMARY STATISTICS ==========')
print('Model parameters:')
print(f'  Wave number k = {waveNumb}')
print(f'  Path loss exponent β = {betaPath}')
print(f'  Absorption coefficient κ = {kappaAbsorp}')
print(f'  Wall distances: a = {a}, b = {b}')
print('\nDeterministic power curve:')
print(f'  Min power = {np.min(sigDetPower):.6f}')
print(f'  Max power = {np.max(sigDetPower):.6f}')
print(f'  Mean power = {np.mean(sigDetPower):.6f}')

if boolePlotOptima and choiceVaryXY < 2:
    print(f'\nTurning points (N = {len(optimaLocPower)}):')
    for i, val in enumerate(optimaLocPower, 1):
        print(f'  {i:2d}: {val:.6f}')

print('\nRandom location model:')
print(f'  Number of samples = {numbRand}')
print(f'  Mean power = {meanSigLocPower:.6f}')
print(f'  Std power = {np.std(sigLocPower):.6f}')
print(f'  Min power = {np.min(sigLocPower):.6f}')
print(f'  Max power = {np.max(sigLocPower):.6f}')
print('========================================\n')

plt.show()
