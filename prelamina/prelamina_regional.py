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

# Map the coordinates to create the curved top profile: y = 1 - 0.4 * exp(-x^2)
coords = mesh.coordinates()
coords[:, 1] *= (1.0 - 0.4 * np.exp(-coords[:, 0]**2))

# ======================================================
# 2. Boundary Markers and Spaces
# ======================================================
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
boundaries.set_all(0)

class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        # Clamp only the bottom wall
        return on_boundary and near(x[1], 0.0)

class Top(SubDomain):
    def inside(self, x, on_boundary):
        # Mark the top curve by excluding the vertical side walls (at x_min and x_max)
        return on_boundary and x[1] > 0.1 and not (near(x[0], x_min) or near(x[0], x_max))

Bottom().mark(boundaries, 13)  # clamp_id
Top().mark(boundaries, 14)    # pressure_id

V = VectorFunctionSpace(mesh, "CG", 2)  #2 mean quadratic

# ======================================================
# 3. Boundary conditions and material parameters
# ======================================================
# Consistent Unit System (Standard for Biomechanics):
# Length: mm | Force: N | Stress/Pressure: MPa (N/mm^2)
# 1 mmHg = 0.0001333 MPa

MPa_TO_mmHg = 1.0 / 0.0001333 # Conversion factor

# Target IOP (will be scaled by load_factor)
iop_target_MPa = 40 * 0.0001333  # Example: 40 mmHg converted to MPa


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
growth_steps = 1000
dt_growth = 0.1
growth_rate = 0.05

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
current_iop_mmHg = iop_target_MPa * MPa_TO_mmHg

  # Define DG0 spaces for local evaluation to prevent memory leaks
V_DG0_tensor = TensorFunctionSpace(mesh, "DG", 0)
V_DG0_scalar = FunctionSpace(mesh, "DG", 0)
atrophy_field = Function(V_DG0_scalar)

for gstep in range(1, growth_steps + 1):
    # --- 1. Element-wise Stress Extraction ---
    F_current = I + grad(u)
    Fe = F_current * inv(Fg)
    
    # Project kinematic tensor
    Fe_proj = project(Fe, V_DG0_tensor, form_compiler_parameters={"quadrature_degree": 2})
    
    # Compute approximate Cauchy stress algebraically
    sigma_approx = mu * (Fe_proj + Fe_proj.T - 2*I) + lmbda * tr(Fe_proj - I)*I
    
    # Project hydrostatic tissue pressure
    p_tissue_expr = -(1.0/3.0) * tr(sigma_approx)
    p_tissue_func = project(p_tissue_expr, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2})
    
    # --- 2. Fast NumPy Biological Math ---
    p_tissue_array = p_tissue_func.vector().get_local()
    
    p_cap_MPa = 25.0 * 0.0001333
    Pt_array = p_cap_MPa - p_tissue_array
    
    P_crit = 5.0 * 0.0001333
    k = 2.0 / 0.0001333
    
    # Calculate local blood flow reduction
    Q_local_array = 0.5 + 0.5 * (1.0 / (1.0 + np.exp(-k * (Pt_array - P_crit))))
    
    # Compute local atrophy trigger (only shrink if flow < 0.8)
    atrophy_array = np.where(Q_local_array < 0.8, 0.8 - Q_local_array, 0.0)
    
    # --- 3. Update FEniCS Field ---
    atrophy_field.vector().set_local(atrophy_array)
    atrophy_field.vector().apply("insert")

    # --- Atrophy/Shrinkage increment ---
    # This updates Fg iteratively using the localized array
    Fg.assign(project((1.0 - growth_rate * atrophy_field * dt_growth) * Fg, Fg_space, form_compiler_parameters={"quadrature_degree": 2}))

    # --- Solve ---
    # The solver 're-reads' the updated Fg values inside the symbolic F_res_growth form
    solve(F_res_growth == 0, u, bcs=[bc_clamp], J=J_res_growth, solver_parameters=solver_params)

    # --- Save ---
    if gstep % 10 == 0:
        avg_atrophy = np.mean(atrophy_array)
        print(f"[Atrophy Step {gstep}] Avg Atrophy Trigger = {avg_atrophy:.6e}, Current IOP = {current_iop_mmHg:.2f} mmHg")

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
