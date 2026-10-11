"""Short science-story captions for the CFD and heat-transfer examples."""

from dataclasses import dataclass


@dataclass(frozen=True)
class Beat:
    """One compact caption in an educational reel."""

    label: str
    title: str
    body: str


@dataclass(frozen=True)
class Story:
    """Narrative and display metadata for one example scene."""

    key: str
    title: str
    hook: str
    units: str
    accent: tuple[int, int, int]
    beats: tuple[Beat, ...]
    takeaway: str
    model_note: str


STORIES: dict[str, Story] = {
    "fixed_end_heat_transfer_1d": Story(
        key="fixed_end_heat_transfer_1d",
        title="Where Does the Heat Go?",
        hook="A hot middle, cooler ends, and one steady heat budget.",
        units="°C · W/m²",
        accent=(70, 220, 235),
        beats=(
            Beat(
                "The setup",
                "Two fixed ends",
                "The rod ends are held at 100°C and 200°C.",
            ),
            Beat(
                "Add a source",
                "Heat from within",
                "Uniform generation adds 100,000 W/m³ throughout the rod.",
            ),
            Beat(
                "Find the field",
                "Steady conduction",
                "The solver balances conductive flux and internal generation.",
            ),
            Beat(
                "Read the colors",
                "Temperature along x",
                "Color marks the computed temperature; rod thickness is illustrative.",
            ),
            Beat(
                "Locate the peak",
                "The warmest point moves",
                "Generation pushes the maximum above either boundary temperature.",
            ),
            Beat(
                "Close the budget",
                "Flux out equals heat in",
                "Both ends together carry away the generated heat per unit area.",
            ),
        ),
        takeaway="At steady state, generated heat leaves through the two fixed-temperature ends.",
        model_note="One-dimensional steady conduction with constant k and uniform volumetric generation.",
    ),
    "convective_end_heat_transfer_1d": Story(
        key="convective_end_heat_transfer_1d",
        title="Two Fluids, One Rod",
        hook="The rod ends are not the same temperature as the fluids around them.",
        units="°C · W/m²K",
        accent=(255, 125, 95),
        beats=(
            Beat(
                "Meet the fluids",
                "Different surroundings",
                "The left and right fluids sit at 100°C and 300°C.",
            ),
            Beat(
                "Exchange heat",
                "Robin conditions",
                "Each end exchanges heat through its own convection coefficient.",
            ),
            Beat(
                "Solve the rod",
                "A steady temperature field",
                "Conduction through the solid links both boundary exchanges.",
            ),
            Beat(
                "Temperatures",
                "Fluid is not solid",
                "The rod-end temperatures differ from the surrounding-fluid temperatures.",
            ),
            Beat(
                "Follow the arrow",
                "Heat flows toward the cold side",
                "The arrow follows the computed direction of net heat flow.",
            ),
            Beat(
                "Check the profile",
                "Compare with analysis",
                "For constant conductivity, the numerical profile matches the analytic solution.",
            ),
        ),
        takeaway="A convection boundary couples a surface temperature to a surrounding fluid; they need not match.",
        model_note="Steady one-dimensional conduction, constant k, with Robin conditions at both ends.",
    ),
    "enhanced_heat_transfer_solver_1d": Story(
        key="enhanced_heat_transfer_solver_1d",
        title="Conductivity That Changes",
        hook="If k depends on temperature, the heat-flow path must adjust too.",
        units="°C · W/mK",
        accent=(180, 135, 255),
        beats=(
            Beat(
                "A new material law",
                "Let k vary with T",
                "The example assumes k(T) = k₀(1 + βT).",
            ),
            Beat(
                "Keep ends convective",
                "Two surrounding fluids",
                "Each end exchanges heat with a fluid at a prescribed temperature.",
            ),
            Beat(
                "Iterate the field",
                "Temperature updates conductivity",
                "The solver repeats until the nonlinear discrete balance converges.",
            ),
            Beat(
                "Track the flux",
                "One heat rate through the rod",
                "Harmonic face conductivity makes neighboring control volumes share flux.",
            ),
            Beat(
                "Read the view",
                "Color still means temperature",
                "The rod shows T; the separate conductivity curve is derived from k(T).",
            ),
            Beat(
                "Check the assumption",
                "The law must stay positive",
                "This chosen linear law is illustrative, not a universal material model.",
            ),
        ),
        takeaway="Temperature-dependent conductivity changes the profile while steady heat flux remains conserved.",
        model_note="Steady conduction without generation; user-selected linear k(T), positive over the solution.",
    ),
    "enhanced_heat_transfer_solver_2d": Story(
        key="enhanced_heat_transfer_solver_2d",
        title="A Heated Square",
        hook="Two insulated sides, two cooled sides, and heat made inside.",
        units="°C · m",
        accent=(255, 165, 75),
        beats=(
            Beat(
                "Set the square",
                "A two-dimensional solid",
                "The field is solved on a rectangular x-y grid.",
            ),
            Beat(
                "Insulate two edges",
                "No heat crosses here",
                "The left and bottom boundaries are insulated.",
            ),
            Beat(
                "Fix the other edges",
                "Cold boundaries draw heat out",
                "The top and right edges are held at prescribed temperatures.",
            ),
            Beat(
                "Generate heat",
                "A source fills the solid",
                "Uniform volumetric generation warms the interior.",
            ),
            Beat(
                "Read the surface",
                "Color is temperature",
                "The raised surface is normalized for display, not physical height.",
            ),
            Beat(
                "Follow the balance",
                "Source meets outward flux",
                "At steady state, generated heat leaves through the fixed-temperature edges.",
            ),
        ),
        takeaway="Insulation redirects heat: generated energy exits through the boundaries that can exchange it.",
        model_note="Steady two-dimensional conduction with uniform generation, insulated left/bottom, fixed top/right.",
    ),
    "heat_convection_solver_1d": Story(
        key="heat_convection_solver_1d",
        title="Advection Meets Diffusion",
        hook="Can a numerical scheme follow a sharp change without smearing it?",
        units="scalar · m/s",
        accent=(80, 210, 160),
        beats=(
            Beat(
                "Set the transport",
                "Flow plus diffusion",
                "A constant velocity carries a scalar while diffusion smooths it.",
            ),
            Beat(
                "Fix both ends",
                "Two prescribed values",
                "The one-dimensional steady problem has Dirichlet endpoint values.",
            ),
            Beat(
                "Choose a scheme",
                "Upwind or central",
                "The discrete convection term changes with the selected scheme.",
            ),
            Beat(
                "Compare curves",
                "Numerical versus analytic",
                "The reference is the steady exponential solution for constant properties.",
            ),
            Beat(
                "Read the rod",
                "Color maps the scalar",
                "The rod field and legend use the numerical and analytic profile range.",
            ),
            Beat(
                "Particle limits",
                "Not a transient solution",
                "Animated random-walk particles illustrate transport; the solved field is steady.",
            ),
        ),
        takeaway="Advection sharpens direction; first-order upwind is robust but adds numerical diffusion.",
        model_note="Steady one-dimensional constant-property advection-diffusion with fixed endpoint values.",
    ),
    "heat_convection_solver_2d": Story(
        key="heat_convection_solver_2d",
        title="Scalar Transport in 2D",
        hook="A prescribed flow carries a scalar across four fixed edges.",
        units="scalar · m/s",
        accent=(85, 190, 255),
        beats=(
            Beat(
                "Choose the field",
                "A rectangular domain",
                "The steady scalar is stored on an x-y grid.",
            ),
            Beat(
                "Prescribe the flow",
                "Velocity is an input",
                "Constant horizontal and vertical components set the transport direction.",
            ),
            Beat(
                "Add diffusion",
                "Smooth across streamlines",
                "Positive diffusion competes with advection in the discrete equation.",
            ),
            Beat(
                "Set all four edges",
                "Fixed scalar values",
                "Dirichlet values are applied around the rectangular boundary.",
            ),
            Beat(
                "Field and arrows",
                "Computed field, prescribed flow",
                "Colors show the scalar; arrows show the input velocity above the surface.",
            ),
            Beat(
                "Model scope",
                "What the model predicts",
                "This solves scalar transport for a prescribed velocity, without computing momentum or pressure.",
            ),
        ),
        takeaway="The steady field reflects the competition between prescribed advection, diffusion, and edge values.",
        model_note="Steady two-dimensional constant-property advection-diffusion with fixed boundary values and upwinding.",
    ),
    "heat_conduction_solver_2d": Story(
        key="heat_conduction_solver_2d",
        title="Heating a Moving Channel",
        hook="Warm fluid downstream, even though the velocity is prescribed.",
        units="°C · W/m²",
        accent=(255, 125, 95),
        beats=(
            Beat(
                "Set the channel",
                "Parallel plates",
                "A prescribed fully developed profile carries fluid through the channel.",
            ),
            Beat(
                "Fix the inlet",
                "Start at one temperature",
                "The inlet temperature is prescribed across the channel height.",
            ),
            Beat(
                "Heat both walls",
                "Equal imposed fluxes",
                "Positive wall heat flux adds energy into the fluid.",
            ),
            Beat(
                "Solve temperature",
                "Development along x",
                "The steady energy equation includes cross-stream and axial conduction.",
            ),
            Beat(
                "Measure the bulk",
                "Weight by velocity",
                "Bulk temperature uses the local velocity profile, not a uniform average.",
            ),
            Beat(
                "Read Nu carefully",
                "A diagnostic, not a flow solve",
                "Nu uses the imposed wall flux and hydraulic diameter 2h; velocity stays prescribed.",
            ),
        ),
        takeaway="The channel computes thermal development in a given flow; it does not solve momentum or pressure.",
        model_note="Steady constant-property energy transport, fixed inlet T, equal wall heat fluxes, zero conductive outlet flux.",
    ),
    "fluid_flow_simulator": Story(
        key="fluid_flow_simulator",
        title="A Lid Sets the Motion",
        hook="Move the top wall, and the cavity flow must respond.",
        units="m/s · kinematic pressure",
        accent=(185, 145, 255),
        beats=(
            Beat(
                "Start with a box",
                "A closed cavity",
                "No-slip walls surround the two-dimensional domain.",
            ),
            Beat(
                "Move the lid",
                "One wall drives flow",
                "The top boundary moves while the other walls remain stationary.",
            ),
            Beat(
                "Iterate momentum",
                "Viscous motion develops",
                "Explicit artificial-compressibility steps update velocity and pressure.",
            ),
            Beat(
                "Read velocity",
                "Magnitude and components",
                "The plots show computed velocity fields and streamlines.",
            ),
            Beat(
                "Trace the circulation",
                "A vortex fills the cavity",
                "Streamlines trace the steady velocity; moving particles are illustrative tracers.",
            ),
            Beat(
                "Know the limits",
                "A teaching approximation",
                "Collocated centered operators can admit pressure checkerboarding.",
            ),
        ),
        takeaway="Boundary motion drives the cavity, while viscosity spreads its influence through the fluid.",
        model_note="Steady approximation using explicit artificial compressibility; pressure may show collocated-grid artifacts.",
    ),
    "obstacle_flow_simulation": Story(
        key="obstacle_flow_simulation",
        title="Flow Around a Block",
        hook="Watch computed flow move around a solid obstacle, step by step.",
        units="lattice velocity · step",
        accent=(255, 155, 85),
        beats=(
            Beat(
                "Set the lattice",
                "A D2Q9 fluid model",
                "Nine discrete velocity directions carry populations across the grid.",
            ),
            Beat(
                "Drive the inlet",
                "Velocity in lattice units",
                "An open x inlet supplies flow; the outlet uses a simple extrapolation.",
            ),
            Beat(
                "Mark the solid",
                "A rectangular obstacle",
                "Halfway bounce-back reflects populations at obstacle links.",
            ),
            Beat(
                "Advance time steps",
                "BGK collision and streaming",
                "Each update relaxes populations then streams them to neighboring nodes.",
            ),
            Beat(
                "Read the snapshots",
                "Velocity, not pressure",
                "Color shows speed; illustrative tracers follow current velocity. Solid cells are masked.",
            ),
            Beat(
                "Periodic edges",
                "Periodic top and bottom",
                "This finite run is not evidence of validated shedding or physical-time behavior.",
            ),
        ),
        takeaway="The lattice method updates discrete populations; its apparent wake depends on grid and boundaries.",
        model_note="Transient D2Q9 BGK with open x boundaries, periodic y boundaries, and a rectangular solid mask.",
    ),
    "transient_heat_diffusion_2d": Story(
        key="transient_heat_diffusion_2d",
        title="Watch Heat Spread",
        hook="Two hot spots on a steel plate. One simple update rule, repeated a thousand times over.",
        units="°C · s",
        accent=(255, 120, 80),
        beats=(
            Beat(
                "Start hot",
                "Two hot spots",
                "A 10 cm steel plate starts near 200°C at two spots; its edges stay at 20°C.",
            ),
            Beat(
                "One small step",
                "Forward in time",
                "Each node's new temperature uses only its four neighbours from the last step.",
            ),
            Beat(
                "Respect the limit",
                "The stability rule",
                "Explicit steps must keep the Fourier number at or below 1/4, or errors explode.",
            ),
            Beat(
                "Spread and merge",
                "Peaks flatten out",
                "Diffusion moves heat from hot to cold, so the two spots widen and join.",
            ),
            Beat(
                "Drain the edges",
                "Cold boundaries win",
                "Heat leaves through the 20°C edges, and the whole plate cools toward them.",
            ),
            Beat(
                "Check the math",
                "Exact decay test",
                "A pure sine mode decays exactly as theory predicts, which verifies the code.",
            ),
        ),
        takeaway="Explicit time stepping is simple and exact enough to verify, but the time step is capped by the grid.",
        model_note="Explicit FTCS solution of the 2D heat equation for carbon steel with fixed-temperature edges.",
    ),
    "rotating_blob_advection_2d": Story(
        key="rotating_blob_advection_2d",
        title="The Blob That Should Not Fade",
        hook="Spin a blob of dye once around. The math says it returns unchanged. The code disagrees.",
        units="scalar · s",
        accent=(186, 140, 255),
        beats=(
            Beat(
                "Set it spinning",
                "Solid-body rotation",
                "The fluid turns like a turntable; one full revolution takes one second.",
            ),
            Beat(
                "Pure transport",
                "No diffusion in the physics",
                "The equation only carries the dye; exactly solved, the peak stays at 1.",
            ),
            Beat(
                "Pick a scheme",
                "First-order upwind",
                "Each face takes its value from the upstream side, which keeps the solution stable.",
            ),
            Beat(
                "Mind the CFL",
                "One cell per step, at most",
                "The time step keeps dye from crossing more than one cell per update.",
            ),
            Beat(
                "Watch it smear",
                "Numerical diffusion",
                "Upwind behaves like added diffusion, so the blob widens and its peak drops.",
            ),
            Beat(
                "Refine the grid",
                "Error shrinks with cells",
                "Halving the cell size roughly halves the smearing of this first-order scheme.",
            ),
        ),
        takeaway="Dye is conserved, but its peak is not: the smearing comes from the scheme, not the physics.",
        model_note="First-order upwind solution of 2D linear advection in solid-body rotation, compared with the exact rotation.",
    ),
}
