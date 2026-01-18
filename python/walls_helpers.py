#!/usr/bin/env python3
"""
walls_helpers.py

Helper functions for electromagnetic signal propagation between parallel walls.

This module provides functions for:
  - Generating transmitter locations (deterministic and random)
  - Computing NLOS signals via method of images
  - Computing NLOS signals via Lerch transcendent function
  - Analyzing signal power and statistics
  - Estimating probability density functions

Reference: "Reflected wireless signals under random spatial sampling"

Author: H. Paul Keeler (with assistance from Claude)
Date: January 2025
"""

import numpy as np
from typing import Tuple, Optional
from scipy.ndimage import uniform_filter1d  # For moving average


# =========================================================================
# TRANSMITTER LOCATION FUNCTIONS
# =========================================================================

def funPosTXPlot(xTX: np.ndarray, yTX: np.ndarray,
                 xDeltaTX: np.ndarray, yDeltaTX: np.ndarray,
                 choiceVaryXY: int, numbPlot: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate deterministic grid of transmitter locations for plotting.

    Creates a grid of evenly-spaced transmitter locations for plotting
    deterministic power curves and surfaces.

    Args:
        xTX: Vector of x-coordinates for transmitter centers [x1, x2, x3]
        yTX: Vector of y-coordinates for transmitter centers [y1, y2, y3]
        xDeltaTX: Vector of x plotting widths [Δx1, Δx2, Δx3]
        yDeltaTX: Vector of y plotting widths [Δy1, Δy2, Δy3]
        choiceVaryXY: Plotting mode (0: vary x; 1: vary y; 2: vary x and y)
        numbPlot: Number of points to generate

    Returns:
        xPlot: X-coordinates for plotting (array)
        yPlot: Y-coordinates for plotting (array)

    Example:
        >>> xTX = np.array([0.25, 0, 0.25])
        >>> yTX = np.array([0, 0.25, 0.25])
        >>> xDeltaTX = np.array([0.2, 0, 0.2])
        >>> yDeltaTX = np.array([0, 0.2, 0.2])
        >>> xPlot, yPlot = funPosTXPlot(xTX, yTX, xDeltaTX, yDeltaTX, 0, 500)
    """
    x = xTX[choiceVaryXY]
    y = yTX[choiceVaryXY]
    xDelta = xDeltaTX[choiceVaryXY]
    yDelta = yDeltaTX[choiceVaryXY]

    if choiceVaryXY == 0:
        # Vary x only, keep y fixed (1D plot)
        xMin = x - xDelta/2
        xMax = x + xDelta/2
        xPlot = np.linspace(xMin, xMax, numbPlot)
        yPlot = y * np.ones(numbPlot)

    elif choiceVaryXY == 1:
        # Keep x fixed, vary y only (1D plot)
        yMin = y - yDelta/2
        yMax = y + yDelta/2
        xPlot = x * np.ones(numbPlot)
        yPlot = np.linspace(yMin, yMax, numbPlot)

    elif choiceVaryXY == 2:
        # Vary both x and y (2D surface plot)
        xMin = x - xDelta/2
        xMax = x + xDelta/2
        yMin = y - yDelta/2
        yMax = y + yDelta/2

        numbMeshX = round(10 * np.sqrt(numbPlot))
        numbMeshY = round(10 * np.sqrt(numbPlot))
        xMesh = np.linspace(xMin, xMax, numbMeshX)
        yMesh = np.linspace(yMin, yMax, numbMeshY)
        xPlot, yPlot = np.meshgrid(xMesh, yMesh)

    else:
        raise ValueError('choiceVaryXY must be 0, 1, or 2')

    return xPlot, yPlot


def funPosTXRand(xTX: np.ndarray, yTX: np.ndarray,
                 xDeltaTX: np.ndarray, yDeltaTX: np.ndarray,
                 choiceVaryXY: int, numbRand: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate random transmitter locations for empirical sampling.

    Randomly samples transmitter locations uniformly within specified regions.
    Used to generate empirical probability density functions (PDFs) from
    Monte Carlo simulation.

    Args:
        xTX: Vector of x-coordinates for transmitter centers [x1, x2, x3]
        yTX: Vector of y-coordinates for transmitter centers [y1, y2, y3]
        xDeltaTX: Vector of x sampling widths [Δx1, Δx2, Δx3]
        yDeltaTX: Vector of y sampling widths [Δy1, Δy2, Δy3]
        choiceVaryXY: Sampling mode (0: vary x; 1: vary y; 2: vary x and y)
        numbRand: Number of random samples to generate

    Returns:
        xRand: Random x-coordinates (array)
        yRand: Random y-coordinates (array)

    Example:
        >>> xRand, yRand = funPosTXRand(xTX, yTX, xDeltaTX, yDeltaTX, 0, 10000)
    """
    x = xTX[choiceVaryXY]
    y = yTX[choiceVaryXY]
    xDelta = xDeltaTX[choiceVaryXY]
    yDelta = yDeltaTX[choiceVaryXY]

    if choiceVaryXY == 0:
        # Vary x only, keep y fixed
        xMin = x - xDelta/2
        xRand = xDelta * np.random.rand(numbRand) + xMin
        yRand = y * np.ones(numbRand)

    elif choiceVaryXY == 1:
        # Keep x fixed, vary y only
        yMin = y - yDelta/2
        xRand = x * np.ones(numbRand)
        yRand = yDelta * np.random.rand(numbRand) + yMin

    elif choiceVaryXY == 2:
        # Vary both x and y (2D rectangular region)
        xMin = x - xDelta/2
        yMin = y - yDelta/2
        xRand = xDelta * np.random.rand(numbRand) + xMin
        yRand = yDelta * np.random.rand(numbRand) + yMin

    else:
        raise ValueError('choiceVaryXY must be 0, 1, or 2')

    return xRand, yRand


# =========================================================================
# NLOS SIGNAL COMPUTATION - DIRECT SUMMATION
# =========================================================================

def funSumNLOS(a: float, b: float, x: np.ndarray, y: np.ndarray,
               betaPath: float, waveNumb: float, kappaAbsorp: float,
               indexFirst: int, indexLast: int) -> np.ndarray:
    """
    Compute non-line-of-sight (NLOS) signal using method of images.

    Calculates the sum of reflected electromagnetic rays between two parallel
    walls using the method of images. This is the direct summation approach.

    Args:
        a: Distance from origin to right wall
        b: Distance from origin to left wall
        x: Transmitter x-coordinates (can be array)
        y: Transmitter y-coordinates (can be array)
        betaPath: Path loss exponent β (signal decays as r^(-β/2))
        waveNumb: Wave number k = 2π/λ
        kappaAbsorp: Absorption coefficient κ ∈ [0,1]
        indexFirst: First reflection order to include (typically 1)
        indexLast: Last reflection order to include (e.g., 100)

    Returns:
        nlosLoc: NLOS signal for random location model (complex array)

    Notes:
        - Uses method of images with 4-case pattern (right-even, left-even,
          left-odd, right-odd)
        - Early termination when terms < 1e-8
        - Each reflection applies absorption factor (-√κ)^n

    Reference:
        Manuscript "Signals reflecting off intelligent walls", Section V
    """
    # Reshape into column vectors
    x = np.atleast_1d(x).flatten()
    y = np.atleast_1d(y).flatten()

    numbSim = len(x)

    # Distance functions for method of images
    funRadRE = lambda ii, x0, y0: np.sqrt(y0**2 + (2*(ii-1)*(a+b) + 2*a - x0)**2)
    funRadLE = lambda ii, x0, y0: np.sqrt(y0**2 + (2*(ii-1)*(a+b) + 2*b + x0)**2)
    funRadRO = lambda kk, x0, y0: np.sqrt(y0**2 + (2*kk*(a+b) + x0)**2)
    funRadLO = lambda kk, x0, y0: np.sqrt(y0**2 + (2*kk*(a+b) - x0)**2)

    # Summation term function
    funSumTerm = lambda ii, kappa0, phase0, rad0: \
        ((-np.sign(kappa0)*np.sqrt(np.abs(kappa0)))**ii) * \
        rad0**(-betaPath/2) * np.exp(1j*phase0)

    # Initialize
    nlosLoc = np.zeros(numbSim, dtype=complex)

    # Determine pattern (cycles through 4 cases)
    numbRad = np.mod(np.arange(indexLast), 4) + 1

    ii_n = 1
    jj_n = 1

    for kk in range(indexFirst, indexLast + 1):
        # Select distance function based on reflection pattern
        if numbRad[kk-1] == 1:  # Right-even
            funRadAny = lambda ii, x0, y0: funRadRE(ii, x0, y0)
        elif numbRad[kk-1] == 2:  # Left-even
            funRadAny = lambda ii, x0, y0: funRadLE(ii, x0, y0)
        elif numbRad[kk-1] == 3:  # Left-odd
            funRadAny = lambda ii, x0, y0: funRadLO(ii, x0, y0)
        elif numbRad[kk-1] == 4:  # Right-odd
            funRadAny = lambda ii, x0, y0: funRadRO(ii, x0, y0)

        # Random location model
        radLoc = funRadAny(ii_n, x, y)
        phaseLoc = waveNumb * radLoc
        termLoc = funSumTerm(jj_n, kappaAbsorp, phaseLoc, radLoc)
        nlosLoc = nlosLoc + termLoc

        # Early termination check
        termAbsMax = np.max(np.abs(termLoc))
        tolTerm = 1e-8

        if termAbsMax < tolTerm:
            break  # Remaining terms will be even smaller

        # Update counters
        ii_n = ii_n + (kk % 4 == 0)
        jj_n = jj_n + (kk % 2 == 0)

    return nlosLoc


# =========================================================================
# NLOS SIGNAL COMPUTATION - LERCH FUNCTION
# =========================================================================

def funTruncLerch(z: complex, s: float, v: np.ndarray) -> np.ndarray:
    """
    Compute the Lerch transcendent function (modified version).

    Calculates a truncated approximation of the modified Lerch transcendent:

        Φ₁(z, s, v) = Σ_{n=1}^∞ z^n / (v + n)^s

    This starts from n=1 instead of n=0, making it appropriate for NLOS
    calculations where the n=0 term (LOS signal) is handled separately.

    Args:
        z: Complex parameter (typically involves absorption and wave number)
        s: Power parameter (typically β/2)
        v: Offset parameter (typically ±r/d), can be scalar or array

    Returns:
        y: Truncated Lerch function value(s), same shape as v

    Notes:
        - Starts from n=1 (excludes LOS term)
        - Uses 400 terms for convergence
        - For standard Lerch Φ (n=0 start), use explicit subtraction

    Reference:
        NIST DLMF: https://dlmf.nist.gov/25.14
    """
    numbFirst = 1      # Start from n=1 (excludes LOS term)
    numbLast = 400     # Sufficient terms for convergence

    nn = np.arange(numbFirst, numbLast + 1)

    v = np.atleast_1d(v)

    if len(v) > 1:
        # Vectorized calculation for multiple v values
        nnMatrix = nn[np.newaxis, :]  # Shape: (1, 400)
        vMatrix = v[:, np.newaxis]    # Shape: (len(v), 1)

        y = np.sum(z**nnMatrix / (vMatrix + nnMatrix)**s, axis=1)

        return y
    else:
        # Scalar calculation
        y = np.sum(z**nn / (v[0] + nn)**s)
        return y


# =========================================================================
# SIGNAL ANALYSIS FUNCTIONS
# =========================================================================

def funSignalPower(sigInput: np.ndarray,
                   expPower: float = 2) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Extract power, real, and imaginary parts of signal.

    Computes the power (magnitude raised to an exponent) and separates the
    real and imaginary components of a complex electromagnetic signal.

    Args:
        sigInput: Complex signal vector or matrix
        expPower: Exponent for power calculation (default: 2)
                  Power = |signal|^expPower

    Returns:
        sigPow: Signal power |sigInput|^expPower
        sigReal: Real part of signal
        sigImag: Imaginary part of signal

    Example:
        >>> S = np.array([1+2j, 3-1j, 2+0j])
        >>> P, Re_S, Im_S = funSignalPower(S, 2)
        >>> # Result: P = [5, 10, 4], Re_S = [1, 3, 2], Im_S = [2, -1, 0]
    """
    sigReal = np.real(sigInput)
    sigImag = np.imag(sigInput)
    sigPow = np.abs(sigInput)**expPower

    return sigPow, sigReal, sigImag


def funSignalStats(sigInput: np.ndarray) -> Tuple[float, float, float, float]:
    """
    Compute mean and variance of complex signal components.

    Calculates basic statistics (mean and variance) separately for the real
    and imaginary parts of a complex electromagnetic signal.

    Args:
        sigInput: Complex signal vector or matrix

    Returns:
        sigMeanReal: Mean of real part
        sigMeanImag: Mean of imaginary part
        sigVarReal: Variance of real part
        sigVarImag: Variance of imaginary part

    Example:
        >>> S = np.random.randn(1000) + 1j*np.random.randn(1000)
        >>> meanRe, meanIm, varRe, varIm = funSignalStats(S)
    """
    sigReal = np.real(sigInput)
    sigImag = np.imag(sigInput)

    sigMeanReal = np.mean(sigReal)
    sigMeanImag = np.mean(sigImag)

    sigVarReal = np.var(sigReal)
    sigVarImag = np.var(sigImag)

    return sigMeanReal, sigMeanImag, sigVarReal, sigVarImag


# =========================================================================
# PDF ESTIMATION FUNCTION
# =========================================================================

def funPDF(X: np.ndarray, booleMethod: bool = True) -> Tuple[np.ndarray, np.ndarray]:
    """
    Estimate probability density function from data.

    Estimates the PDF of data using either histogram binning with smoothing
    or kernel density estimation.

    Args:
        X: Data vector from which to estimate PDF
        booleMethod: PDF estimation method
                    True: Histogram with 100 bins + 3-point moving average
                    False: Kernel density estimation (requires scipy.stats)

    Returns:
        f: Estimated probability density values
        xi: Corresponding x-values where density is evaluated

    Notes:
        - Method True: Faster but may miss fine details
        - Method False: Smoother but requires scipy.stats.gaussian_kde
        - Histogram method recommended for large datasets (N > 10000)

    Example:
        >>> X = np.random.randn(10000)
        >>> f, xi = funPDF(X, True)
        >>> plt.plot(xi, f)
    """
    if booleMethod:
        # METHOD 1: Histogram with smoothing
        pdfNEmp, vecX = np.histogram(X, bins=100, density=True)

        # Apply 3-point moving average for smoothing
        pdfNEmp = uniform_filter1d(pdfNEmp, size=3, mode='nearest')

        # Output
        f = pdfNEmp
        xi = vecX[:-1]  # Bin left edges

    else:
        # METHOD 0: Kernel density estimation
        from scipy.stats import gaussian_kde
        kde = gaussian_kde(X)
        xi = np.linspace(np.min(X), np.max(X), 100)
        f = kde(xi)

    return f, xi
