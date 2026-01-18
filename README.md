# WallsFading

Code repository for the manuscript:

**"Reflected wireless signals under random spatial sampling"**  
H. Paul Keeler  


## Summary

This repository provides MATLAB and Python implementations demonstrating that **randomly positioning a transmitter between parallel reflecting walls produces singularities in the probability density function (PDF)** of received signal power. These singularities occur at power values corresponding to turning points in the deterministic power function—a fundamental phenomenon arising from the coupling between wave interference and spatial geometry.

## Key Finding

When a transmitter is uniformly distributed in a region bounded by reflecting walls:

1. The deterministic signal power P(r) oscillates due to constructive/destructive interference
2. At turning points where P'(r) = 0, the standard change-of-variables formula fails
3. This failure manifests as **inverse square-root singularities** in the PDF of received power
4. Multiple spatial turning points can collapse into single PDF singularities due to many-to-one mapping

This contrasts sharply with the classical random phase model, which produces smooth (often Rayleigh-like) distributions.

## Repository Structure

```
WallsFading/
├── README.md                 # This file
├── LICENSE
├── matlab/
│   ├── README.md             # MATLAB-specific documentation
│   ├── walls_fading_figures.m    # Main script: generates Figures 8-12
│   ├── RandomRays.m              # Interactive exploration tool
│   ├── CompareLerchVsSummation.m # Validation script
│   ├── funSumNLOS.m              # Method of images (direct summation)
│   ├── funTruncLerch.m           # Lerch transcendent function
│   ├── funPosTXPlot.m            # Deterministic position grid
│   ├── funPosTXRand.m            # Random position sampling
│   ├── funSignalPower.m          # Power extraction
│   ├── funSignalStats.m          # Signal statistics
│   └── funPDF.m                  # PDF estimation
└── python/
    ├── README.md             # Python-specific documentation
    ├── walls_fading_figures.py   # Main script: generates Figures 8-12
    ├── RandomRays.py             # Interactive exploration tool
    ├── CompareLerchVsSummation.py # Validation script
    └── walls_helpers.py          # Helper functions module
```

## Quick Start

### MATLAB

```matlab
cd matlab
run('walls_fading_figures.m')   % Generate Figures 8-12
```

### Python

```bash
cd python
pip install numpy scipy matplotlib
python walls_fading_figures.py   # Generate Figures 8-12
```

## Parameters

Default parameters used throughout (matching manuscript Figures 5-14):

| Parameter | Symbol | Value | Description |
|-----------|--------|-------|-------------|
| Wall separation | d | 1.0 | Distance between walls (a = b = 0.5) |
| Path loss exponent | β | 4.0 | Signal attenuation exponent |
| Reflection coefficient | κ | 0.5 | Power fraction reflected per bounce |
| Monte Carlo samples | N | 50,000 | Random locations for PDF estimation |

Wave numbers vary by figure:

| Figure | k | Configuration |
|--------|-----|---------------|
| 8 | 100 | Vary x ∈ [0.15, 0.35], y = 0 |
| 9 | 1000 | Vary y ∈ [0.15, 0.35], x = 0 |
| 10 | 100 | Vary y ∈ [-0.1, 0.6], x = 0 |
| 11 | 100 | Vary y ∈ [0.1, 0.6], x = 0 |
| 12 | 10 | Vary x ∈ [0.05, 0.45], y ∈ [-0.2, 0.2] |

## Mathematical Background

### Signal Model

For a transmitter at distance r from the origin between two parallel walls separated by distance d, the NLOS signal is expressed using the Lerch transcendent (Equation 21 in manuscript):

```
S(r) = (e^{-jkr} / d^{β/2}) [Φ(ζ, β/2, -r/d) - correction]
     + (e^{+jkr} / d^{β/2}) [Φ(ζ, β/2, +r/d) - correction]
```

where:
- Φ(ζ, s, α) = Σ_{n=0}^∞ ζⁿ / (n + α)^s is the Lerch transcendent
- ζ = -√κ · e^{jkd} encodes absorption and round-trip phase
- The correction terms subtract the n=0 (line-of-sight) contribution

### PDF Singularity Mechanism

At a turning point t where the power function satisfies P'(t) = 0, the PDF exhibits (Proposition IV.3):

```
f_V(v) ~ C / √|v - P(t)|   as v → P(t)
```

This inverse square-root singularity arises because the Jacobian of the transformation vanishes at turning points, causing probability to concentrate near these power values.

### Key Decomposition

The signal can be decomposed as S(r) = A(r)·W(r), where:
- A(r) represents slowly-varying attenuation
- W(r) captures rapid wave interference oscillations

The singularities arise from turning points in |W(r)|², demonstrating that the phenomenon is fundamentally geometric rather than dependent on specific propagation details.

## Validation

Both implementations have been cross-validated:

- **Lerch vs Direct Summation**: Agreement to ~10⁻¹³ relative error
- **MATLAB vs Python**: Matching results for all figures
- **Convergence**: 200-400 terms sufficient for k ≤ 1000

Run the validation scripts:
```matlab
% MATLAB
run('CompareLerchVsSummation.m')
```
```bash
# Python
python CompareLerchVsSummation.py
```

## Output

Each figure contains two panels:
- **Top**: Deterministic power curve P(x,y) with turning points marked (red dots)
- **Bottom**: Empirical PDF from Monte Carlo sampling with singularity locations (red dashed lines)

## Requirements

### MATLAB
- MATLAB R2017b or later
- Signal Processing Toolbox (for `islocalmin`, `islocalmax`)
- Image Processing Toolbox (optional, for 2D extrema in Figure 12)

### Python
- Python 3.7+
- NumPy
- SciPy
- Matplotlib

## Applications

The theoretical results have implications for:
- **Intelligent surface deployment**: Understanding statistical behavior of signals in environments with reconfigurable intelligent surfaces (RIS)
- **Network optimization**: Predicting coverage probability in indoor/canyon environments
- **Stochastic geometry**: Bridging deterministic ray-tracing with statistical fading analysis

## References

1. NIST Digital Library of Mathematical Functions, §25.14 (Lerch Transcendent)
2. Lindell & Alanen, "Exact image theory for the Sommerfeld half-space problem" (method of images)
3. Björnson & Sanguinetti, "Rayleigh fading modeling and channel hardening for reconfigurable intelligent surfaces"

## License

See [LICENSE](LICENSE) file.

## Citation

If you use this code in your research, please cite:

```bibtex
@article{keeler2026wallsfading,
  title={Reflected wireless signals under random spatial sampling},
  author={Keeler, H. Paul},
  journal={Print},
  year={2026},
  note={Submitted}
}
```

