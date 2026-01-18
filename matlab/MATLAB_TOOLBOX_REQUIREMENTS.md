# MATLAB Toolbox Requirements

## Overview

This document details the MATLAB toolbox dependencies for the "Reflected wireless signals" code.

## Required Toolboxes

### Signal Processing Toolbox ✓ REQUIRED

**Used by:**
- `walls_fading_figures.m` (main figure generation script)
- `RandomRays.m` (interactive demonstration)

**Functions used:**
- `islocalmin()` - Find local minima in power curves
- `islocalmax()` - Find local maxima in power curves

**Purpose:** Identifying turning points (local extrema) in the deterministic power function P(x,y), which correspond to singularities in the empirical PDF.

**Alternative:** If you don't have this toolbox, you could implement manual peak detection:

```matlab
% Manual replacement for islocalmin/islocalmax
function idx = manual_localmin(data)
    n = length(data);
    idx = false(n, 1);
    for i = 2:n-1
        if data(i) < data(i-1) && data(i) < data(i+1)
            idx(i) = true;
        end
    end
end
```

However, the toolbox version is more robust and handles edge cases better.

---

### Image Processing Toolbox (Optional)

**Used by:**
- `RandomRays.m` (only for 2D mode with `choiceVaryXY = 3`)

**Functions used:**
- `imregionalmin()` - Find regional minima in 2D power surfaces
- `imregionalmax()` - Find regional maxima in 2D power surfaces

**Purpose:** Identifying turning points in 2D surfaces when both x and y vary simultaneously.

**Workaround:** If you don't have this toolbox:
1. Only use 1D modes in `RandomRays.m` (set `choiceVaryXY = 1` or `2`)
2. The main figure generation script `walls_fading_figures.m` doesn't use 2D, so it will work without this toolbox

**Note:** The manuscript figures (8-12) are all 1D or can be generated without 2D turning point detection, so this toolbox is truly optional.

---

## Base MATLAB Functions (No Toolbox Required)

The following functions are part of base MATLAB and require no additional toolboxes:

- `linspace()` - Generate linearly spaced vectors
- `rand()`, `randn()` - Random number generation
- `rng()` - Random number generator control
- `histcounts()` - Histogram computation
- `histogram()` - Histogram plotting
- `plot()`, `subplot()`, `surf()` - Plotting functions
- `exp()`, `abs()`, `sqrt()` - Mathematical functions
- `meshgrid()` - Create 2D grids
- `hypot()` - Euclidean distance
- All standard arithmetic and logical operations

---

## Checking Your Installation

To check which toolboxes you have installed:

```matlab
ver
```

Look for:
```
Signal Processing Toolbox                      Version X.X
Image Processing Toolbox                       Version X.X
```

Or check programmatically:

```matlab
% Check Signal Processing Toolbox
if license('test', 'Signal_Toolbox')
    fprintf('✓ Signal Processing Toolbox: Available\n');
else
    fprintf('✗ Signal Processing Toolbox: NOT available\n');
    fprintf('  → Required for turning point detection\n');
end

% Check Image Processing Toolbox
if license('test', 'Image_Toolbox')
    fprintf('✓ Image Processing Toolbox: Available\n');
else
    fprintf('✗ Image Processing Toolbox: NOT available\n');
    fprintf('  → Only needed for 2D mode in RandomRays.m\n');
end
```

---

## Minimum Configuration

**To run main figure generation (`walls_fading_figures.m`):**
- Base MATLAB R2018b or later
- Signal Processing Toolbox

**To run validation (`CompareLerchVsSummation.m`):**
- Base MATLAB only (no toolboxes required!)

**To run interactive demo (`RandomRays.m`):**
- Signal Processing Toolbox (for 1D modes)
- Image Processing Toolbox (additionally, for 2D mode)

---

## Installation Instructions

If you're missing required toolboxes:

### Option 1: Via MATLAB Add-On Explorer (Recommended)
1. Open MATLAB
2. Click "Home" tab → "Add-Ons" → "Get Add-Ons"
3. Search for "Signal Processing Toolbox"
4. Click "Install" (requires license)

### Option 2: Via MathWorks Account
1. Log in to your MathWorks account
2. Go to "My Account" → "My Software"
3. Select your MATLAB license
4. Check the boxes for required toolboxes
5. Download and install

### Option 3: Contact Your Administrator
If using an institutional license, contact your IT department or MATLAB administrator to request toolbox access.

---

## Cost Considerations

### Academic/Student Licenses
Most universities provide MATLAB with all toolboxes at no additional cost. Check with your institution's IT department.

### Individual Licenses
- Signal Processing Toolbox: ~$1000 (perpetual) or included in many bundles
- Image Processing Toolbox: ~$1000 (perpetual) or included in many bundles
- Alternatively, consider the Python version which uses only free, open-source packages

---

## Function Mapping: MATLAB vs Python

If toolbox costs are prohibitive, consider using the Python implementation instead:

| MATLAB Function | Toolbox | Python Equivalent | Package |
|----------------|---------|-------------------|---------|
| `islocalmin()` | Signal Proc. | `scipy.signal.argrelextrema(..., np.less)` | scipy (free) |
| `islocalmax()` | Signal Proc. | `scipy.signal.argrelextrema(..., np.greater)` | scipy (free) |
| `imregionalmin()` | Image Proc. | Manual implementation | numpy (free) |
| `imregionalmax()` | Image Proc. | Manual implementation | numpy (free) |

The Python implementation provides identical results with zero licensing costs.

---

## Summary Table

| Component | Toolbox Required | Used By | Can Work Without? |
|-----------|------------------|---------|-------------------|
| Main figures (8-12) | Signal Processing | walls_fading_figures.m | No |
| Validation script | None | CompareLerchVsSummation.m | Yes - no toolboxes needed |
| Interactive demo (1D) | Signal Processing | RandomRays.m | No |
| Interactive demo (2D) | Signal + Image | RandomRays.m | Use 1D mode instead |
| Helper functions | None | fun*.m | Yes - all base MATLAB |

---

## Recommendations

1. **If you have Signal Processing Toolbox:** You're good to go! All main features work.

2. **If you only have base MATLAB:** 
   - Run validation script (no toolboxes needed)
   - Use Python implementation for figures
   - Or implement manual peak detection

3. **If you have neither toolbox:**
   - Use the Python implementation (free, identical results)
   - Or request toolbox access from your institution

4. **For new users:** 
   - Start with Python (free, easier to install)
   - Move to MATLAB if integration with existing MATLAB workflows is needed

---

## Questions?

If you have questions about toolbox requirements or alternatives, please contact:

H. Paul Keeler  
Email: hpkeeler@unimelb.edu.au
