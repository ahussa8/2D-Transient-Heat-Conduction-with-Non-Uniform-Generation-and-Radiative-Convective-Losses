# 2D Heat Conduction with Radiative–Convective Surface Exchange

This project models transient 2D heat conduction in a diamond plate with internal volumetric heat generation, coupled with surface radiative–convective heat loss analysis. Implemented using MATLAB.

## 🔧 Methodology

### Part I: Unsteady Heat Conduction
- Solved the 2D heat conduction PDE using implicit finite difference with Gauss–Seidel iteration.
- Included two internal heat sources (15,000 and 10,000 W/m³) on a 10×10 grid.
- Time-varying Dirichlet BCs on top/left/right and insulated Neumann BC at the bottom.
- Displayed transient temperature fields at 0, 0.7, and 1.4s.

### Part II: Radiative–Convective Exchange
- Used surface temperature-dependent emissivity to model radiosity and q_rad.
- Computed q_conv using flat-plate Nusselt number correlations.
- Evaluated net surface flux (q_net) across 300–1000 K.
- Plotted blackbody and graybody spectral radiance using Planck’s law.

## 📊 Results

- Peak temperatures reached 360 K and 345 K near internal sources after 1.4s.
- q_net crossed zero near 800 K, indicating thermal equilibrium.
- Radiative losses dominated above 1200 K.
- Spectral radiance plots confirmed expected trends from Wien's law.
