#!/usr/bin/env python3
"""
CompareLerchVsSummation.py

Demonstrates mathematical equivalence between two methods for computing
NLOS (non-line-of-sight) electromagnetic signals between parallel walls:
  1. Direct summation using method of images
  2. Closed-form expression using modified Lerch transcendent function

The modified Lerch function Φ₁(ζ,s,α) is defined as:

  Φ₁(ζ,s,α) = Σ_{n=1}^∞ ζ^n / (n + α)^s

This differs from the standard Lerch Φ(ζ,s,α) by starting at n=1 instead
of n=0. The n=0 term represents the direct line-of-sight (LOS) signal,
which is handled separately, making Φ₁ appropriate for NLOS calculations.

REFERENCE:
  Manuscript "Reflected wireless signals under random spatial sampling"
  Proposition V.1 (Equation 29): Closed-form NLOS signal expression
"""

import numpy as np
import matplotlib.pyplot as plt
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
from walls_helpers import funSumNLOS, funTruncLerch

print('='*40)
print('Comparing Lerch Function vs Direct Sum')
print('='*40)
print()

# =========================================================================
# PARAMETERS (Using manuscript Figure 9 settings)
# =========================================================================

# Physical parameters
a = 0.5                 # Distance to right wall
b = 0.5                 # Distance to left wall
d = a + b               # Total wall separation
betaPath = 4.0          # Path loss exponent β
kappaAbsorp = 0.5       # Absorption coefficient κ
waveNumb = 100          # Wave number k

# Summation parameters
indexFirst = 1          # Start from first reflection (n=1)
indexLast = 200         # Use 200 terms for convergence

# Test points
nTest = 10              # Number of test points
xTest = np.linspace(0.15, 0.35, nTest)  # x coordinates
yTest = np.zeros(nTest)                  # y = 0 (fixed)

# =========================================================================
# METHOD 1: Direct Summation (Method of Images)
# =========================================================================

print('METHOD 1: Direct summation using method of images')
print('-' * 50)

startTime = time.time()
nlosDirect = funSumNLOS(a, b, xTest, yTest, betaPath,
                         waveNumb, kappaAbsorp, indexFirst, indexLast)
timeDirect = time.time() - startTime

print(f'Computed {nTest} test points in {timeDirect:.4f} seconds')
print(f'Average per point: {timeDirect/nTest:.6f} seconds\n')

# =========================================================================
# METHOD 2: Closed-Form Using Modified Lerch Function Φ₁
# =========================================================================

print('METHOD 2: Closed-form using modified Lerch function Φ₁')
print('-' * 50)

# Preallocate
nlosLerch = np.zeros(nTest, dtype=complex)

startTime = time.time()
for i in range(nTest):
    # Distance from origin
    r = np.hypot(xTest[i], yTest[i])
    if r < 0.001:
        r = 0.001  # Avoid singularity

    # Common parameters for Lerch function
    zeta = -np.sqrt(kappaAbsorp) * np.exp(1j * waveNumb * d)
    s = betaPath / 2

    # First term: right images (Equation 29a in manuscript)
    alpha1 = -r / d
    phi1 = funTruncLerch(zeta, s, alpha1)
    term1 = (np.exp(-1j * waveNumb * r) / (d**s)) * phi1

    # Second term: left images (Equation 29b in manuscript)
    alpha2 = r / d
    phi2 = funTruncLerch(zeta, s, alpha2)
    term2 = (np.exp(1j * waveNumb * r) / (d**s)) * phi2

    # Total NLOS signal
    nlosLerch[i] = term1 + term2

timeLerch = time.time() - startTime

print(f'Computed {nTest} test points in {timeLerch:.4f} seconds')
print(f'Average per point: {timeLerch/nTest:.6f} seconds')
print(f'Speedup factor: {timeDirect/timeLerch:.2f}x faster\n')

# =========================================================================
# COMPARISON: Compute Differences
# =========================================================================

print('COMPARISON RESULTS')
print('-' * 50)

# Compute differences
diffComplex = nlosDirect - nlosLerch
diffMagnitude = np.abs(diffComplex)
diffPower = np.abs(nlosDirect)**2 - np.abs(nlosLerch)**2

# Statistical summary
print('Complex signal differences:')
print(f'  Max |S_direct - S_Lerch| = {np.max(diffMagnitude):.2e}')
print(f'  Mean |S_direct - S_Lerch| = {np.mean(diffMagnitude):.2e}')
print(f'  Relative error = {100*np.max(diffMagnitude)/np.mean(np.abs(nlosDirect)):.2e}%\n')

print('Power differences:')
print(f'  Max |P_direct - P_Lerch| = {np.max(np.abs(diffPower)):.2e}')
print(f'  Mean |P_direct - P_Lerch| = {np.mean(np.abs(diffPower)):.2e}')
print(f'  Relative error = {100*np.max(np.abs(diffPower))/np.mean(np.abs(nlosDirect)**2):.2e}%\n')

# =========================================================================
# DETAILED COMPARISON TABLE
# =========================================================================

print('DETAILED POINT-BY-POINT COMPARISON')
print('-' * 50)
print('Point |   x   |   Direct Method   |   Lerch Method    | Difference')
print('------|-------|-------------------|-------------------|-----------')

for i in range(nTest):
    r = np.hypot(xTest[i], yTest[i])
    magDirect = np.abs(nlosDirect[i])
    magLerch = np.abs(nlosLerch[i])
    diff = np.abs(nlosDirect[i] - nlosLerch[i])

    print(f'{i+1:5d} | {xTest[i]:.3f} | '
          f'{np.real(nlosDirect[i]):8.6f} + {np.imag(nlosDirect[i]):8.6f}i | '
          f'{np.real(nlosLerch[i]):8.6f} + {np.imag(nlosLerch[i]):8.6f}i | '
          f'{diff:.2e}')
print()

# =========================================================================
# VISUALIZATION: Signal Magnitude and Power
# =========================================================================

fig = plt.figure(figsize=(12, 8))
fig.suptitle('Lerch vs Direct Summation Comparison', fontsize=16, fontweight='bold')

# Compute signal magnitudes and power
magDirect = np.abs(nlosDirect)
magLerch = np.abs(nlosLerch)
powerDirect = magDirect**2
powerLerch = magLerch**2

# Subplot 1: Signal magnitude comparison
ax1 = plt.subplot(2, 2, 1)
plt.plot(xTest, magDirect, 'bo-', linewidth=2, markersize=8, label='Direct Sum')
plt.plot(xTest, magLerch, 'rx--', linewidth=2, markersize=10, label='Lerch Φ₁')
plt.xlabel('x position', fontsize=16)
plt.ylabel('Signal magnitude |S|', fontsize=16)
plt.title('Signal Magnitude: Two Methods', fontsize=16)
plt.legend(loc='best', fontsize=12)
plt.grid(True, alpha=0.3)

# Subplot 2: Magnitude difference (zoomed)
ax2 = plt.subplot(2, 2, 2)
plt.plot(xTest, diffMagnitude, 'ko-', linewidth=1.5, markersize=6)
plt.xlabel('x position', fontsize=16)
plt.ylabel('|S_direct - S_Lerch|', fontsize=16)
plt.title('Magnitude Difference (max error shown)', fontsize=16)
plt.grid(True, alpha=0.3)
plt.ylim([0, np.max(diffMagnitude)*1.2])
plt.text(np.mean(xTest), np.max(diffMagnitude)*0.9,
         f'Max error: {np.max(diffMagnitude):.2e}',
         ha='center', fontsize=12, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

# Subplot 3: Power comparison
ax3 = plt.subplot(2, 2, 3)
plt.plot(xTest, powerDirect, 'bo-', linewidth=2, markersize=8, label='Direct Sum')
plt.plot(xTest, powerLerch, 'rx--', linewidth=2, markersize=10, label='Lerch Φ₁')
plt.xlabel('x position', fontsize=16)
plt.ylabel('Signal power P = |S|²', fontsize=16)
plt.title('Signal Power: Two Methods', fontsize=16)
plt.legend(loc='best', fontsize=12)
plt.grid(True, alpha=0.3)

# Subplot 4: Power difference
ax4 = plt.subplot(2, 2, 4)
plt.plot(xTest, np.abs(diffPower), 'ko-', linewidth=1.5, markersize=6)
plt.xlabel('x position', fontsize=16)
plt.ylabel('|P_direct - P_Lerch|', fontsize=16)
plt.title('Power Difference', fontsize=16)
plt.grid(True, alpha=0.3)
plt.ylim([0, np.max(np.abs(diffPower))*1.2])
plt.text(np.mean(xTest), np.max(np.abs(diffPower))*0.9,
         f'Max error: {np.max(np.abs(diffPower)):.2e}',
         ha='center', fontsize=12, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.tight_layout()

# =========================================================================
# MATHEMATICAL EXPLANATION
# =========================================================================

print('MATHEMATICAL EQUIVALENCE')
print('='*40)
print()

print('The NLOS signal can be expressed as an infinite sum over reflections:\n')
print('  S_NLOS(r) = Σ_{n=1}^∞ (-√κ)^n × ℓ(r̂_n) × e^(jkr̂_n)')
print('              + Σ_{n=1}^∞ (-√κ)^n × ℓ(t̂_n) × e^(jkt̂_n)\n')

print('where r̂_n and t̂_n are distances to right and left image sources,')
print('and ℓ(r) = r^(-β/2) is the path loss function.\n')

print('For the symmetric case (a = b = d/2), this simplifies to:\n')
print('  S_NLOS(r) = (e^(-jkr) / d^(β/2)) × [Φ₁(ζ, β/2, -r/d)]')
print('              + (e^(+jkr) / d^(β/2)) × [Φ₁(ζ, β/2, +r/d)]\n')

print('where:')
print('  ζ = -√κ × e^(jkd)  (absorption and phase per round-trip)')
print('  Φ₁(ζ,s,α) = Σ_{n=1}^∞ ζ^n/(n+α)^s  (modified Lerch function)\n')

print('KEY INSIGHT:')
print('-' * 12)
print('The modified Lerch Φ₁ starts from n=1 (not n=0) because:')
print('  • The n=0 term represents the direct LOS path')
print('  • NLOS calculations exclude this direct path')
print('  • This is mathematically equivalent to using standard Lerch Φ')
print('    with explicit subtraction: Φ(ζ,s,α) - α^(-s)\n')

print(f'Convergence: Both methods use {indexLast} terms. Increasing to 300-400 terms')
print('may improve accuracy for large wave numbers (k > 1000).\n')

# =========================================================================
# CONCLUSION
# =========================================================================

print('CONCLUSION')
print('='*40)
print()

if np.max(diffMagnitude) < 1e-6:
    print('✓ VALIDATION SUCCESSFUL\n')
    print('The two methods produce identical results within numerical precision:')
    print(f'  • Maximum difference: {np.max(diffMagnitude):.2e}')
    print(f'  • Relative error: {100*np.max(diffMagnitude)/np.mean(magDirect):.2e}%\n')
    print('The closed-form Lerch expression (Proposition V.1) is confirmed to be')
    print('mathematically equivalent to the direct method-of-images summation.')
else:
    print('⚠ WARNING: Significant difference detected\n')
    print('Consider increasing number of terms (indexLast) for better convergence.')

print()
print('='*40)

plt.show()
