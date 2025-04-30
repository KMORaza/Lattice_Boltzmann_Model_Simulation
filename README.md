## Simulation of Lattice-Boltzmann Model written in JavaFX

* The application simulates fluid dynamics using the Lattice-Boltzmann model and helps to visualize fluid flow with different lattice types, collision operators, and boundary conditions.
* Supported lattice types :-
  * D2Q9: 2D lattice with 9 velocity directions, 100x100 grid.
  * D3Q15: 3D lattice with 15 velocity directions, 50x50x50 grid.
  * D3Q19: 3D lattice with 19 velocity directions, 50x50x50 grid.
  * D3Q27: 3D lattice with 27 velocity directions, 50x50x50 grid.
* Each lattice defines velocity vectors (C) and weights (W) for particle distributions. The simulation computes fluid density (rho) and velocities (ux, uy, uz) at each grid point, stored in arrays.
* Functions performed by simulation loop :-
  * Compute macroscopic quantities: Calculates density and velocities from the distribution functions (f).
  * Collide: Updates particle distributions based on the chosen collision operator.
  * Stream: Moves distributions along velocity directions.
  * Apply boundary conditions: Enforces rules at the grid boundaries.
  * Visualization: Updates the canvas with the current velocity field.
* Collision operators :-
  * SRT (Single Relaxation Time): Uses a single relaxation parameter (omega = 1/0.6) for all moments, applied to all lattices.
  * MRT (Multiple Relaxation Time): Uses multiple relaxation rates for different moments, implemented for D2Q9 only.
  * TRT (Two Relaxation Time): Uses symmetric and antisymmetric relaxation times, for D2Q9 only.
  * Entropic: Adjusts relaxation to maximize entropy, for D2Q9 only. For 3D lattices, non-SRT operators default to SRT.
* Boundary conditions :-
  * Bounce-back: Reflects particles at walls, simulating no-slip conditions.
  * Velocity: Sets a fixed velocity (0.1 in x-direction) at the left boundary (x=0), with bounce-back elsewhere.
  * Pressure: Fixes density (1.0) at the right boundary (x=NX-1), with bounce-back elsewhere.
  * Periodic: Wraps particles around opposite boundaries.
  * Inlet/Outlet: Combines velocity at the left and pressure at the right, with bounce-back on top/bottom.
  * Open: Copies distributions from neighboring cells, allowing free flow through boundaries.
* The simulation initializes with a uniform density (1.0) and zero velocity, except for specific boundary conditions (e.g., 0.1 x-velocity at x=0 for Velocity/Inlet/Outlet). An obstacle is placed to demonstrate flow interaction, visible in the visualization.

---

| ![](https://github.com/KMORaza/Lattice_Boltzmann_Model_Simulation/blob/main/src/screenshots/screen01.png) | ![](https://github.com/KMORaza/Lattice_Boltzmann_Model_Simulation/blob/main/src/screenshots/screen02.png) |
|-----------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------|
| ![](https://github.com/KMORaza/Lattice_Boltzmann_Model_Simulation/blob/main/src/screenshots/screen03.png) | ![](https://github.com/KMORaza/Lattice_Boltzmann_Model_Simulation/blob/main/src/screenshots/screen04.png) |
| ![](https://github.com/KMORaza/Lattice_Boltzmann_Model_Simulation/blob/main/src/screenshots/screen05.png) | ![](https://github.com/KMORaza/Lattice_Boltzmann_Model_Simulation/blob/main/src/screenshots/screen06.png) |
| ![](https://github.com/KMORaza/Lattice_Boltzmann_Model_Simulation/blob/main/src/screenshots/screen07.png) | ![](https://github.com/KMORaza/Lattice_Boltzmann_Model_Simulation/blob/main/src/screenshots/screen08.png) |
