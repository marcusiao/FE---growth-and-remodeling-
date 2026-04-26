from dolfin import *
import numpy as np
import os

# ======================================================
# 1. Mesh Generation (2D Prelaminar Tissue)
# ======================================================
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Define domain resolution and extents
nx, ny = 60, 40
x_min, x_max = -2.0, 2.0
y_min, y_max = 0.0, 1.0

# Create a base rectangle mesh
mesh = RectangleMesh(Point(x_min, y_min), Point(x_max, y_max), nx, ny)

# Map the coordinates to create the curved top profile: y = 1 - 0.2 * exp(-x^2)
coords = mesh.coordinates()
coords[:, 1] *= (1.0 - 0.2 * np.exp(-coords[:, 0]**2))

# ======================================================
# 2. Boundary Markers and Spaces
# ======================================================
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
boundaries.set_all(0)

class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        # The bottom edge remains at y=0 after the transformation
        return on_boundary and near(x[1], 0.0)

class Top(SubDomain):
    def inside(self, x, on_boundary):
        # Mark the top curve by excluding the vertical side walls (at x_min and x_max)
        return on_boundary and x[1] > 0.1 and not (near(x[0], x_min) or near(x[0], x_max))

Bottom().mark(boundaries, 13) # clamp_id
Top().mark(boundaries, 14)    # pressure_id

V = VectorFunctionSpace(mesh, "CG", 2)

# ======================================================
# 3. Boundary conditions and material parameters
# ======================================================
# Consistent Unit System (Standard for Biomechanics):
# Length: mm | Force: N | Stress/Pressure: MPa (N/mm^2)
# 1 mmHg = 0.0001333 MPa

# Blood flow atrophy parameters
MAP_mmHg = 90.0 # Assuming constant Mean Arterial Pressure
HEALTHY_IOP_mmHg = 10.0 # Baseline healthy IOP
OPP_BASELINE_mmHg = (2/3) * MAP_mmHg - HEALTHY_IOP_mmHg # Calculated baseline healthy OPP
OPP_THRESHOLD_REDUCTION_mmHg = 20.0 # OPP reduction up to this value causes no blood flow change
OPP_MAX_REDUCTION_FOR_50_PERCENT_FLOW_mmHg = 35.0 # OPP reduction at which blood flow drops to 50%
BLOOD_FLOW_DROP_AT_MAX_OPP_REDUCTION = 0.5 # 50% drop in blood flow
MPa_TO_mmHg = 1.0 / 0.0001333 # Conversion factor

# Target IOP (will be scaled by load_factor)
iop_target_MPa = 25 * 0.0001333  # Example: 40 mmHg converted to MPa


zero = Constant((0.0, 0.0))
clamp_id = 13
bc_clamp = DirichletBC(V, zero, boundaries, clamp_id)

nu = 0.499                                  # Poisson's ratio (near incompressible)
mu_val = 62.3 / 1000                        # Float value for Python calculations, convert kPa to MPa
E_val = 2 * mu_val * (1 + nu)               # Back cal the Young's mod
lmbda_val = (E_val * nu) / ((1 + nu) * (1 - 2 * nu)) # Calculate lmbda as a float

mu = Constant(mu_val)                       # Wrapped for FEniCS efficiency
lmbda = Constant(lmbda_val)                 # Wrapped for FEniCS efficiency

u = Function(V)
du = TrialFunction(V)
v = TestFunction(V)

d = mesh.geometry().dim()
I = Identity(d)

# Define a custom measure with fixed quadrature degree for consistency
dx_c = Measure("dx", domain=mesh, metadata={"quadrature_degree": 4})

# ======================================================
# 4. Loading Phase (step load)
# ======================================================
load_steps = 100
pressure_id = 14 # ID for IOP boundary
load_factor = Constant(0.0)
ds_measure = Measure("ds", domain=mesh, subdomain_data=boundaries)
n = FacetNormal(mesh)

F = I + grad(u)
C = F.T * F
Ic = tr(C)
J = det(F)
psi = (mu / 2) * (Ic - 3) - mu * ln(J) + (lmbda / 2) * (ln(J))**2   #neo hookean equation

# Potential energy with normal pressure: Pi = Psi - (-P * n . u) = Psi + P * n . u
Pi_total = psi * dx_c + load_factor * iop_target_MPa * dot(n, u) * ds_measure(pressure_id)
F_res = derivative(Pi_total, u, v)
J_res = derivative(F_res, u, du)

# Shared solver parameters for robustness
solver_params = {"newton_solver": {"relative_tolerance": 1e-8, 
                                   "absolute_tolerance": 1e-10, 
                                   "maximum_iterations": 25, 
                                   "linear_solver": "mumps"}}

xdmf_load = XDMFFile(mesh.mpi_comm(), "u_loading_phase.xdmf")
xdmf_load.parameters["flush_output"] = True
xdmf_load.parameters["functions_share_mesh"] = True

print("=== Loading Phase ===")
for step in range(1, load_steps + 1):
    load_factor.assign(step / load_steps) # load_factor now scales from 0 to 1
    solve(F_res == 0, u, bcs=[bc_clamp], J=J_res, solver_parameters=solver_params)
    xdmf_load.write(u, float(step))
xdmf_load.close()
print("✅ Loading phase finished and saved.")

# ======================================================
# 5. Remodeling Phase (maintain load)
# ======================================================
growth_steps = 100
dt_growth = 0.1
growth_rate = 0.01

Fg_space = TensorFunctionSpace(mesh, "CG", 1)
Fg = Function(Fg_space)
Fg.assign(project(Identity(d), Fg_space))

# Keep load at final value
load_factor.assign(1.0) # Keep load at full target_load (1.0 * iop_target_MPa)

xdmf_growth = XDMFFile(mesh.mpi_comm(), "u_growth_phase.xdmf")
xdmf_growth.parameters["flush_output"] = True
xdmf_growth.parameters["functions_share_mesh"] = True

xdmf_Fg = XDMFFile(mesh.mpi_comm(), "Fg_growth_phase.xdmf")
xdmf_Fg.parameters["flush_output"] = True
xdmf_Fg.parameters["functions_share_mesh"] = True

# --- Define variational problem once outside the loop ---
F_current_sym = I + grad(u)
Fe_sym = F_current_sym * inv(Fg)  # This is the elastic part Fe from F = Fe * Fg
C_sym = Fe_sym.T * Fe_sym
Ic_sym = tr(C_sym)
J_sym = det(Fe_sym)
psi_growth = (mu/2)*(Ic_sym - 3) - mu*ln(J_sym) + (lmbda/2)*(ln(J_sym))**2 # Use current mu, lmbda
Pi_total_growth = psi_growth*dx_c + load_factor * iop_target_MPa * dot(n, u) * ds_measure(pressure_id)
F_res_growth = derivative(Pi_total_growth, u, v)
J_res_growth = derivative(F_res_growth, u, du)

print("=== Remodeling Phase (Atrophy/Shrinkage) ===")
# The atrophy is triggered by the steady-state IOP reached after loading
current_iop_mmHg = iop_target_MPa * MPa_TO_mmHg

# --- Calculate atrophy trigger value once (it is constant during this phase) ---
OPP_current_mmHg = (2/3) * MAP_mmHg - current_iop_mmHg
OPP_reduction_mmHg = OPP_BASELINE_mmHg - OPP_current_mmHg

blood_flow_reduction = 0.0
if OPP_reduction_mmHg > OPP_THRESHOLD_REDUCTION_mmHg:
    OPP_reduction_beyond_threshold = OPP_reduction_mmHg - OPP_THRESHOLD_REDUCTION_mmHg
    
    # Slope of blood flow reduction (0.5 flow drop over 15 mmHg OPP reduction)
    slope = BLOOD_FLOW_DROP_AT_MAX_OPP_REDUCTION / (OPP_MAX_REDUCTION_FOR_50_PERCENT_FLOW_mmHg - OPP_THRESHOLD_REDUCTION_mmHg)
    
    blood_flow_reduction = slope * OPP_reduction_beyond_threshold
    
    # Cap blood flow reduction at 50%
    blood_flow_reduction = min(blood_flow_reduction, BLOOD_FLOW_DROP_AT_MAX_OPP_REDUCTION)

atrophy_trigger_value = blood_flow_reduction

for gstep in range(1, growth_steps + 1):
    # --- Atrophy/Shrinkage increment ---
    # This updates Fg iteratively: Fg_new = (1 - alpha) * Fg_old
    # Even if atrophy_trigger_value is constant, Fg shrinks exponentially.
    Fg.assign(project((1.0 - growth_rate * atrophy_trigger_value * dt_growth) * Fg, Fg_space))

    # --- Solve ---
    # The solver 're-reads' the updated Fg values inside the symbolic F_res_growth form
    solve(F_res_growth == 0, u, bcs=[bc_clamp], J=J_res_growth, solver_parameters=solver_params)

    # --- Save ---
    if gstep % 10 == 0:
        print(f"[Atrophy Step {gstep}] Atrophy Trigger = {atrophy_trigger_value:.6e}, Current IOP = {current_iop_mmHg:.2f} mmHg, OPP = {OPP_current_mmHg:.2f} mmHg")

    xdmf_growth.write(u, float(gstep))
    if gstep % 20 == 0:
        xdmf_Fg.write(Fg, float(gstep))

xdmf_growth.close()
xdmf_Fg.close()
print("✅ Remodeling phase finished and saved.")

# ======================================================
# 6. Unloading Phase (gradual release from remodeled state)
# ======================================================
unload_steps = 100
xdmf_unload = XDMFFile(mesh.mpi_comm(), "u_unloading_phase.xdmf")
xdmf_unload.parameters["flush_output"] = True
xdmf_unload.parameters["functions_share_mesh"] = True

# --- Define variational problem once for unloading ---
F_current_un = I + grad(u)
Fe_un = F_current_un * inv(Fg)
C_un = Fe_un.T * Fe_un
Ic_un = tr(C_un)
J_un = det(Fe_un)
psi_un = (mu / 2) * (Ic_un - 3) - mu * ln(J_un) + (lmbda / 2) * (ln(J_un))**2
Pi_total_unload = psi_un * dx_c + load_factor * iop_target_MPa * dot(n, u) * ds_measure(pressure_id)
F_res_unload = derivative(Pi_total_unload, u, v)
J_res_unload = derivative(F_res_unload, u, du)

print("=== Unloading Phase ===")
for step in range(1, unload_steps + 1):
    # Gradually reduce load
    load_factor.assign(1.0 * (1 - step / unload_steps)) # Scale from 1.0 down to 0.0

    # --- Solve equilibrium ---
    solve(F_res_unload == 0, u, bcs=[bc_clamp], J=J_res_unload, solver_parameters=solver_params)

    # --- Save displacement ---
    xdmf_unload.write(u, float(step))

xdmf_unload.close()
print("✅ Unloading phase finished and saved.")
