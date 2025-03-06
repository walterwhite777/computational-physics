import numpy as np
import matplotlib.pyplot as plt

# Constants
hbar = 1.0  # Reduced Planck constant (in atomic units)
m = 1.0     # Electron mass (in atomic units)

# Define the potential barrier and well
def repulsive_barrier(x, V0=10, a=1):
    return V0 * (np.abs(x) <= a)

def attractive_well(x, V0=-10, a=1):
    return V0 * (np.abs(x) <= a)

# Define the Schrödinger equation in 1D
def schrodinger_equation(energy, potential, x):
    dx = x[1] - x[0]
    psi = np.zeros_like(x)
    psi[0] = 0.0
    psi[1] = 0.01  # Small initial value
    for i in range(1, len(x) - 1):
        psi[i + 1] = (2 * (1 - (dx ** 2 / hbar**2) * (potential[i] - energy) * m) * psi[i] - psi[i - 1])
    return psi

# Define spatial grid
x = np.linspace(-5, 5, 1000)

# Compute wave functions for scattering states
energy = 5.0  # Energy of the incoming electron
barrier_potential = repulsive_barrier(x)
well_potential = attractive_well(x)

psi_barrier = schrodinger_equation(energy, barrier_potential, x)
psi_well = schrodinger_equation(energy, well_potential, x)

# Normalize wave functions
psi_barrier /= np.sqrt(np.trapz(psi_barrier**2, x))
psi_well /= np.sqrt(np.trapz(psi_well**2, x))

# Plot results
plt.figure(figsize=(10, 5))
plt.plot(x, psi_barrier, label='Scattering by Repulsive Barrier')
plt.plot(x, psi_well, label='Scattering by Attractive Well')
plt.xlabel('x')
plt.ylabel('Wave Function')
plt.title('Electron Scattering States Comparison')
plt.legend()
plt.grid(True)

# Plot potential barriers for comparison
plt.figure(figsize=(10, 5))
plt.plot(x, barrier_potential, label='Repulsive Barrier')
plt.plot(x, well_potential, label='Attractive Well')
plt.xlabel('x')
plt.ylabel('Potential Energy')
plt.title('Potential Energy Comparison')
plt.legend()
plt.grid(True)

plt.show()