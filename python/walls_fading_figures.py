"""
Reproduce Figures 8-12 from manuscript "Reflected wireless signals under random spatial sampling"
Random location model with empirical PDF estimation
Python version for cross-validation with MATLAB implementation
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema

# Set random seed for reproducibility
np.random.seed(42)

print('Generating Figures 8-12 with beta = 4...\n')

# Parameters
d = 1.0          # a = b = 0.5, so d = 1.0
beta = 4.0       # Path loss exponent (CORRECTED from 2.0)
kappa = 0.5      # Reflection coefficient
N_samples = 50000  # Number of Monte Carlo samples
N_plot = 500     # Number of points for plotting

def lerch_transcendent(zeta, s, alpha, max_terms=300):
    """Compute Lerch transcendent Phi(zeta, s, alpha)"""
    y = 0.0 + 0.0j
    for n in range(max_terms):
        denom = n + alpha
        if abs(denom) < 1e-12:
            continue
        term = (zeta**n) / (denom**s)
        y += term
        if abs(term) < 1e-12:
            break
    return y

def nlos_signal_lerch(r, d, kappa, k, beta, max_terms=200):
    """NLOS signal using Lerch function (Equation 21 from manuscript)"""
    if abs(r) < 1e-10:
        r = 1e-10
    
    # First term (right images)
    zeta1 = -np.sqrt(kappa) * np.exp(1j * k * d)
    alpha1 = -r / d
    phi1 = lerch_transcendent(zeta1, beta/2, alpha1, max_terms)
    term1 = (np.exp(-1j * k * r) / (d**(beta/2))) * (phi1 - (-d/r)**(beta/2))
    
    # Second term (left images)
    zeta2 = -np.sqrt(kappa) * np.exp(1j * k * d)
    alpha2 = r / d
    phi2 = lerch_transcendent(zeta2, beta/2, alpha2, max_terms)
    term2 = (np.exp(1j * k * r) / (d**(beta/2))) * (phi2 - (d/r)**(beta/2))
    
    return term1 + term2

def find_turning_points(x_vals, P_vals):
    """Find local minima and maxima"""
    min_idx = argrelextrema(P_vals, np.less)[0]
    max_idx = argrelextrema(P_vals, np.greater)[0]
    turning_idx = np.sort(np.concatenate([min_idx, max_idx]))
    return turning_idx

# ========================================================================
# FIGURE 8: Random location in X direction, k=100
# ========================================================================
print('Generating Figure 8 (k=100, vary x)...')

k = 100
x_min, x_max = 0.15, 0.35

# Top plot: Power function
x_range = np.linspace(x_min, x_max, N_plot)
P_power = np.zeros(N_plot)

for i, x in enumerate(x_range):
    r = max(abs(x), 0.001)
    S = nlos_signal_lerch(r, d, kappa, k, beta)
    P_power[i] = np.abs(S)**2

# Find turning points
turning_idx = find_turning_points(x_range, P_power)
turning_powers = P_power[turning_idx]

# Bottom plot: Monte Carlo sampling
print(f'  Sampling {N_samples} random locations...')
X_samples = np.random.uniform(x_min, x_max, N_samples)
P_samples = np.zeros(N_samples)

for i, x in enumerate(X_samples):
    r = max(abs(x), 0.001)
    S = nlos_signal_lerch(r, d, kappa, k, beta)
    P_samples[i] = np.abs(S)**2
    if (i + 1) % 10000 == 0:
        print(f'    Processed {i+1}/{N_samples}')

# Create figure
fig, axes = plt.subplots(2, 1, figsize=(10, 10))

axes[0].plot(x_range, P_power, 'b-', linewidth=2)
axes[0].plot(x_range[turning_idx], P_power[turning_idx], 'ro', markersize=8, linewidth=2)
axes[0].set_xlabel('x', fontsize=12)
axes[0].set_ylabel('Power value P(x,y)', fontsize=12)
axes[0].set_title('Vary x, keep y fixed.', fontsize=12)
axes[0].grid(True)
axes[0].set_xlim([x_min, x_max])

axes[1].hist(P_samples, bins=100, density=True, color=[0.4, 0.6, 1], edgecolor='none')
for tp in turning_powers:
    axes[1].axvline(tp, color='r', linestyle='--', linewidth=1.5)
axes[1].set_xlabel('Power value s', fontsize=12)
axes[1].set_ylabel('Empirical PDF of P', fontsize=12)
axes[1].grid(True)

plt.tight_layout()
plt.savefig('/home/claude/Figure8_python.png', dpi=150)
print('  Saved Figure8_python.png')
plt.close()

# ========================================================================
# FIGURE 9: Random location in Y direction, k=1000
# ========================================================================
print('\nGenerating Figure 9 (k=1000, vary y)...')

k = 1000
y_min, y_max = 0.15, 0.35

# Top plot: Power function
y_range = np.linspace(y_min, y_max, N_plot)
P_power = np.zeros(N_plot)

for i, y in enumerate(y_range):
    r = max(abs(y), 0.001)
    S = nlos_signal_lerch(r, d, kappa, k, beta, max_terms=300)
    P_power[i] = np.abs(S)**2

turning_idx = find_turning_points(y_range, P_power)
turning_powers = P_power[turning_idx]

# Bottom plot: Monte Carlo sampling
print(f'  Sampling {N_samples} random locations...')
Y_samples = np.random.uniform(y_min, y_max, N_samples)
P_samples = np.zeros(N_samples)

for i, y in enumerate(Y_samples):
    r = max(abs(y), 0.001)
    S = nlos_signal_lerch(r, d, kappa, k, beta, max_terms=300)
    P_samples[i] = np.abs(S)**2
    if (i + 1) % 10000 == 0:
        print(f'    Processed {i+1}/{N_samples}')

fig, axes = plt.subplots(2, 1, figsize=(10, 10))

axes[0].plot(y_range, P_power, 'b-', linewidth=2)
axes[0].plot(y_range[turning_idx], P_power[turning_idx], 'ro', markersize=8, linewidth=2)
axes[0].set_xlabel('y', fontsize=12)
axes[0].set_ylabel('Power value P(x,y)', fontsize=12)
axes[0].set_title('Vary y, keep x fixed.', fontsize=12)
axes[0].grid(True)
axes[0].set_xlim([y_min, y_max])

axes[1].hist(P_samples, bins=100, density=True, color=[0.4, 0.6, 1], edgecolor='none')
for tp in turning_powers:
    axes[1].axvline(tp, color='r', linestyle='--', linewidth=1.5)
axes[1].set_xlabel('Power value s', fontsize=12)
axes[1].set_ylabel('Empirical PDF of P', fontsize=12)
axes[1].grid(True)

plt.tight_layout()
plt.savefig('/home/claude/Figure9_python.png', dpi=150)
print('  Saved Figure9_python.png')
plt.close()

# ========================================================================
# FIGURE 10: Random location in Y, wider range, k=100
# ========================================================================
print('\nGenerating Figure 10 (k=100, wider y range)...')

k = 100
y_min, y_max = -0.5, 0.5

y_range = np.linspace(y_min, y_max, N_plot)
P_power = np.zeros(N_plot)

for i, y in enumerate(y_range):
    r = max(abs(y), 0.001)
    S = nlos_signal_lerch(r, d, kappa, k, beta)
    P_power[i] = np.abs(S)**2

turning_idx = find_turning_points(y_range, P_power)
turning_powers = P_power[turning_idx]

print(f'  Sampling {N_samples} random locations...')
Y_samples = np.random.uniform(y_min, y_max, N_samples)
P_samples = np.zeros(N_samples)

for i, y in enumerate(Y_samples):
    r = max(abs(y), 0.001)
    S = nlos_signal_lerch(r, d, kappa, k, beta)
    P_samples[i] = np.abs(S)**2
    if (i + 1) % 10000 == 0:
        print(f'    Processed {i+1}/{N_samples}')

fig, axes = plt.subplots(2, 1, figsize=(10, 10))

axes[0].plot(y_range, P_power, 'b-', linewidth=2)
axes[0].plot(y_range[turning_idx], P_power[turning_idx], 'ro', markersize=8, linewidth=2)
axes[0].set_xlabel('y', fontsize=12)
axes[0].set_ylabel('Power value P(x,y)', fontsize=12)
axes[0].set_title('Vary y, keep x fixed.', fontsize=12)
axes[0].grid(True)
axes[0].set_xlim([y_min, y_max])

axes[1].hist(P_samples, bins=100, density=True, color=[0.4, 0.6, 1], edgecolor='none')
for tp in turning_powers:
    axes[1].axvline(tp, color='r', linestyle='--', linewidth=1.5)
axes[1].set_xlabel('Power value s', fontsize=12)
axes[1].set_ylabel('Empirical PDF of P', fontsize=12)
axes[1].grid(True)

plt.tight_layout()
plt.savefig('/home/claude/Figure10_python.png', dpi=150)
print('  Saved Figure10_python.png')
plt.close()

# ========================================================================
# FIGURE 11: Random location in Y, smaller interval, k=100
# ========================================================================
print('\nGenerating Figure 11 (k=100, narrower y range)...')

k = 100
y_min, y_max = 0, 0.6

y_range = np.linspace(y_min, y_max, N_plot)
P_power = np.zeros(N_plot)

for i, y in enumerate(y_range):
    r = max(abs(y), 0.001)
    S = nlos_signal_lerch(r, d, kappa, k, beta)
    P_power[i] = np.abs(S)**2

turning_idx = find_turning_points(y_range, P_power)
turning_powers = P_power[turning_idx]

print(f'  Sampling {N_samples} random locations...')
Y_samples = np.random.uniform(y_min, y_max, N_samples)
P_samples = np.zeros(N_samples)

for i, y in enumerate(Y_samples):
    r = max(abs(y), 0.001)
    S = nlos_signal_lerch(r, d, kappa, k, beta)
    P_samples[i] = np.abs(S)**2
    if (i + 1) % 10000 == 0:
        print(f'    Processed {i+1}/{N_samples}')

fig, axes = plt.subplots(2, 1, figsize=(10, 10))

axes[0].plot(y_range, P_power, 'b-', linewidth=2)
axes[0].plot(y_range[turning_idx], P_power[turning_idx], 'ro', markersize=8, linewidth=2)
axes[0].set_xlabel('y', fontsize=12)
axes[0].set_ylabel('Power value P(x,y)', fontsize=12)
axes[0].set_title('Vary y, keep x fixed.', fontsize=12)
axes[0].grid(True)
axes[0].set_xlim([y_min, y_max])

axes[1].hist(P_samples, bins=100, density=True, color=[0.4, 0.6, 1], edgecolor='none')
for tp in turning_powers:
    axes[1].axvline(tp, color='r', linestyle='--', linewidth=1.5)
axes[1].set_xlabel('Power value s', fontsize=12)
axes[1].set_ylabel('Empirical PDF of P', fontsize=12)
axes[1].grid(True)

plt.tight_layout()
plt.savefig('/home/claude/Figure11_python.png', dpi=150)
print('  Saved Figure11_python.png')
plt.close()

# ========================================================================
# FIGURE 12: Random location in both X and Y, k=10
# ========================================================================
print('\nGenerating Figure 12 (k=10, vary x and y)...')

k = 10
x_min, x_max = 0.05, 0.45
y_min, y_max = -0.2, 0.2

# Top plot: 3D surface
N_mesh = 40
x_vals = np.linspace(x_min, x_max, N_mesh)
y_vals = np.linspace(y_min, y_max, N_mesh)
X_grid, Y_grid = np.meshgrid(x_vals, y_vals)
P_surface = np.zeros_like(X_grid)

print('  Computing 3D surface...')
for i in range(X_grid.shape[0]):
    for j in range(X_grid.shape[1]):
        r = np.sqrt(X_grid[i, j]**2 + Y_grid[i, j]**2)
        r = max(r, 0.001)
        S = nlos_signal_lerch(r, d, kappa, k, beta)
        P_surface[i, j] = np.abs(S)**2

# Bottom plot: Random 2D sampling
print(f'  Sampling {N_samples} random 2D locations...')
X_samples = np.random.uniform(x_min, x_max, N_samples)
Y_samples = np.random.uniform(y_min, y_max, N_samples)
P_samples = np.zeros(N_samples)

for i in range(N_samples):
    r = np.sqrt(X_samples[i]**2 + Y_samples[i]**2)
    r = max(r, 0.001)
    S = nlos_signal_lerch(r, d, kappa, k, beta)
    P_samples[i] = np.abs(S)**2
    if (i + 1) % 10000 == 0:
        print(f'    Processed {i+1}/{N_samples}')

fig = plt.figure(figsize=(12, 10))

ax1 = fig.add_subplot(2, 1, 1, projection='3d')
surf = ax1.plot_surface(X_grid, Y_grid, P_surface, cmap='viridis')
ax1.set_xlabel('x', fontsize=11)
ax1.set_ylabel('y', fontsize=11)
ax1.set_zlabel('P(x,y)', fontsize=11)
ax1.set_title('Vary x and y', fontsize=12)
fig.colorbar(surf, ax=ax1, shrink=0.5)

ax2 = fig.add_subplot(2, 1, 2)
ax2.hist(P_samples, bins=100, density=True, color=[0.4, 0.6, 1], edgecolor='none')
ax2.set_xlabel('s', fontsize=12)
ax2.set_ylabel('Empirical PDF of P', fontsize=12)
ax2.set_title('Random location', fontsize=12)
ax2.grid(True)

plt.tight_layout()
plt.savefig('/home/claude/Figure12_python.png', dpi=150)
print('  Saved Figure12_python.png')
plt.close()

print('\n' + '='*60)
print('All figures generated successfully with beta = 4!')
print('='*60)
print('\nGenerated files:')
print('  - Figure8_python.png   (X random, k=100)')
print('  - Figure9_python.png   (Y random, k=1000)')
print('  - Figure10_python.png  (Y random, wide range, k=100)')
print('  - Figure11_python.png  (Y random, narrow range, k=100)')
print('  - Figure12_python.png  (X,Y both random, k=10)')
print('\nParameters used:')
print(f'  d = {d}, beta = {beta}, kappa = {kappa}')
print(f'  N_samples = {N_samples}')
