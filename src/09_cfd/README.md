# CFD and heat-transfer demonstrations

These examples pair small numerical models with Matplotlib and VTK views. The
heat models have analytical and conservation checks. The flow models are teaching
implementations; their boundary rules and discrete equations are tested, but they
have not been validated against a grid-converged external flow benchmark.

Run an example from the repository root:

```bash
python -m src --list cfd
python -m src heat_conduction_solver_2d
```

| Example | Problem actually computed | How to interpret the visualization |
| --- | --- | --- |
| [Transient heat diffusion, 2D](transient_heat_diffusion_2d.py) *(beginner)* | Explicit FTCS time marching of the 2D heat equation in a 10 cm carbon-steel plate with two Gaussian hot spots and fixed 20°C edges. The step respects `Fo_x + Fo_y <= 1/2`. | Height and color are temperature at recorded times; snapshots are denser early, where the field changes fastest. Exceed the step limit to see the instability. |
| [Rotating blob advection, 2D](rotating_blob_advection_2d.py) *(beginner)* | Explicit first-order upwind solution of linear advection in solid-body rotation (one revolution per second) at CFL 0.8, with zero inflow at the edges. | The exact solution rotates the blob unchanged; the lower, wider numerical blob is numerical diffusion. Height and color are concentration. |
| [Fixed ends, 1D](fixed_end_heat_transfer_1d.py) | Steady conduction with constant conductivity, uniform volumetric generation, and two fixed end temperatures. | Rod colors are the computed temperature at each axial position; displayed rod radius is illustrative. |
| [Convective ends, 1D](convective_end_heat_transfer_1d.py) | Steady conduction with two Robin conditions describing exchange with surrounding fluids. | Fluid temperatures differ from rod end temperatures. Heat-flow arrows follow the sign of the boundary heat exchange. |
| [Variable conductivity, 1D](enhanced_heat_transfer_solver_1d.py) | Steady conduction with `k(T)=k_c*(1+beta*T)` and convective ends, without internal generation. | Rod color shows temperature; the conductivity profile is a derived field. Nonpositive conductivity is rejected. |
| [Heated square, 2D](enhanced_heat_transfer_solver_2d.py) | Steady conduction in a square with uniform generation, insulated left/bottom edges, and fixed top/right temperatures. | Color is temperature. The raised surface is normalized temperature, not a three-dimensional physical domain. |
| [Scalar transport, 1D](heat_convection_solver_1d.py) | Steady advection–diffusion with constant velocity/diffusion and fixed endpoint values, comparing discretization schemes with an analytical solution. | Numerical and analytical profiles are distinct. Animated particles illustrate transport; they do not solve the transient scalar equation. |
| [Scalar transport, 2D](heat_convection_solver_2d.py) | Steady advection–diffusion with constant prescribed `(u,v)`, positive diffusion, fixed edge values, and upwind convection. | Color and contours show the computed scalar; height is normalized. Arrows show prescribed velocity on a plane above the surface. Particles are illustrative tracers. |
| [Heated channel, 2D](heat_conduction_solver_2d.py) | Steady energy transport in prescribed laminar parallel-plate flow, mean speed 1 m/s, fixed inlet temperature, equal wall heat fluxes, and zero conductive outlet flux. | This computes temperature, not fluid velocity/pressure. Colors use physical x/y coordinates. Bulk temperature is velocity weighted; `Nu` uses hydraulic diameter `2*h`. |
| [Lid-driven cavity](fluid_flow_simulator.py) | A steady approximation of viscous incompressible flow, using explicit artificial-compressibility iterations. | Velocity plots show computed magnitude/components and streamlines. Kinematic pressure has an arbitrary reference and may contain collocated-grid artifacts. |
| [Obstacle flow](obstacle_flow_simulation.py) | Transient D2Q9 BGK lattice-Boltzmann flow, open in x, periodic in y, with a rectangular solid obstacle. | Frames are post-step velocity snapshots in lattice units. Solid points are masked. The color scale stays fixed across frames. Top/bottom are periodic boundaries, not solid walls. |

## Numerical conventions and checks

Two-dimensional fields use `field[y_index, x_index]`. VTK point ordering is
x-fastest, and each attached scalar belongs to that point's actual coordinates.
Color legends use the same range as their mapper. Decorative display height is
explicitly labeled. It does not change scalar values.

Finite-volume coefficients include each face's length per unit depth. Omitting
these lengths changes the relative x/y diffusion on rectangular grids. Flux
boundaries use half control volumes (quarter at two intersecting insulated
edges), and their source and boundary equations participate in convergence.
This follows the face-flux balance described in the
[NIST FiPy discretization documentation](https://pages.nist.gov/fipy/en/stable/numerical/discret.html).

The 2D heat iterations check the equation residual scaled by the diagonal
coefficient, in temperature/scalar units. The nonlinear 1D model checks a
relative residual of the assembled interior and boundary equations. They have finite iteration limits and raise
`ConvergenceError` on failure. Shared TDMA validates its inputs, avoids mutating
callers, and detects unsafe pivots. The cavity checks both divergence and
momentum residual, rather than accepting small velocity updates as proof of a
steady incompressible solution.

Regression checks include:

- FTCS exact discrete decay of a sine mode, the continuous decay rate on a fine grid, blow-up past the Fourier limit, and the maximum principle.
- Upwind CFL step size, monotonicity, flux-form dye conservation, first-order refinement toward the exact rotated blob, and peak loss from numerical diffusion.
- Fixed-end and convective-end analytical temperature solutions and boundary heat balances.
- Nonlinear conductive flux balance, constant-conductivity reduction, and refinement toward an integrated-conductivity solution.
- Stable analytical 1D transport at large positive and negative Peclet numbers.
- Harmonic and affine fields on rectangular grids; upwind transport refinement toward an exponential solution.
- Second-order refinement of the insulated/fixed square toward a separated cosine-series solution, and generation/outward-flux balance.
- Channel no-slip forward velocity, insulated uniform temperature, heating/cooling signs, and full discrete energy balance.
- Developed equal-flux channel `Nu ≈ 140/17 ≈ 8.2353`. A wall derivative is a separate refinement diagnostic; the heat-transfer coefficient uses the imposed wall flux.
- Lattice equilibrium mass/momentum, non-equilibrium inlet reconstruction, streaming, obstacle reflection, and finite short runs.
- Rectangular field-to-VTK coordinate/scalar mapping, quad orientation, and legend ranges.

Run the checks without opening interactive windows:

```bash
xvfb-run -a python -m unittest discover -s tests -v
```

On a desktop with a display, omit `xvfb-run -a`.

## Limits of the models

In the generic scalar models, `Gamma` is density times scalar diffusivity.
For temperature transport with `Gamma=k`, the parameter called `rho` must be
volumetric heat capacity (density times specific heat). Alternatively, use mass
density for `rho` and `Gamma=k/c_p`. The heated-channel model explicitly includes
`rho*c`, so its `k` is thermal conductivity.

Upwind convection is first order and adds numerical diffusion. A finite residual
does not establish mesh independence; refine the grid before interpreting a
sharp thermal or scalar boundary layer quantitatively. The nonlinear model
assumes its linear conductivity law remains positive over the solution.

Channel velocity is an imposed fully developed profile. Thermal development is
computed, with constant properties, axial conduction, and no buoyancy or
velocity/temperature coupling. Both wall fluxes are positive **into the fluid**
for heating. The inlet fixes temperature even at its wall corners. Wall heat is
integrated over the active control volumes, starting at `x=dx/2`. Heat-transfer
coefficients are undefined where wall and bulk temperatures coincide, so the
inlet coefficient is displayed as missing (`NaN`).

The cavity uses collocated centered pressure/divergence operators without a
Rhie–Chow correction. They can admit checkerboard pressure modes. Its explicit
step is constrained by advection, artificial acoustic propagation, and viscous
diffusion. Passing the discrete residual checks is useful for this model; it
does not make its pressure a validated engineering result.

The obstacle example defaults to `Re=100`: the full 520×180, 3000-step run
remains finite with this setting. The former `Re=220` default became unstable
at the inlet after hundreds of steps with the corrected boundary formulation.
Higher Reynolds numbers remain configurable; invalid populations raise an error
rather than being clipped into a plausible-looking frame.

The obstacle model uses the inlet reconstruction from
[Zou and He](https://arxiv.org/abs/comp-gas/9611001), halfway reflection at solid
links, and a simple extrapolation outlet. Its Reynolds number is based on the
actual obstacle height, `Re=U_in*obstacle_height/nu`. BGK near its relaxation
limit can become unstable; a finite snapshot sequence alone does not establish
correct drag, shedding frequency, or an asymptotic wake. The outlet and periodic
y boundaries influence the result. No physical seconds or meters are inferred
from lattice units.
