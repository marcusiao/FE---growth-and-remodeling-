from dolfin import *
import numpy as np
import meshio
import os

# ======================================================
# 1. Mesh Generation (2D Prelaminar Tissue)
# ======================================================
os.chdir(os.path.dirname(os.path.realpath(__file__)))

file_name = "PT_model"

# --- Convert .msh to .xdmf if not already present ---
if not (os.path.exists(f"{file_name}_domain.xdmf") and os.path.exists(f"{file_name}_boundaries.xdmf")):
    msh = meshio.read(f"{file_name}.msh")
    cells_vol = msh.get_cells_type("triangle")
    cell_data_vol = msh.get_cell_data("gmsh:physical", "triangle")
    meshio.write_points_cells(
        f"{file_name}_domain.xdmf",
        msh.points[:, :2],
        [("triangle", cells_vol)],
        cell_data={"name_to_read": [cell_data_vol]}
    )

    cells_surf = msh.get_cells_type("line")
    cell_data_surf = msh.get_cell_data("gmsh:physical", "line")
    meshio.write_points_cells(
        f"{file_name}_boundaries.xdmf",
        msh.points[:, :2],
        [("line", cells_surf)],
        cell_data={"name_to_read": [cell_data_surf]}
    )

# ======================================================
# 2. Load mesh and define spaces
# ======================================================
mesh = Mesh()
with XDMFFile(f"{file_name}_domain.xdmf") as infile:
    infile.read(mesh)

mvc_bnd = MeshValueCollection("size_t", mesh, mesh.topology().dim() - 1)
with XDMFFile(f"{file_name}_boundaries.xdmf") as infile:
    infile.read(mvc_bnd, "name_to_read")
boundaries = MeshFunction("size_t", mesh, mvc_bnd)

ID_LC_FLOOR = 1
ID_SCLERA_WALL = 2
ID_CUP_SURFACE_IOP = 3
ID_SYMMETRY_AXIS = 4
ID_RETINA_BOUNDARY = 5

V = VectorFunctionSpace(mesh, "CG", 2)  #2 mean quadratic

# ======================================================
# 3. Boundary conditions and material parameters
# ======================================================
# Consistent Unit System (Standard for Biomechanics):
# Length: mm | Force: N | Stress/Pressure: MPa (N/mm^2)
# 1 mmHg = 0.0001333 MPa

# Target IOP and CSFP (will be scaled by load_factor)
iop_target_MPa = 40 * 0.0001333  # Increased from 40 to 150 mmHg to increase tissue stress
csfp_target_MPa = 0.00 * 0.0001333 # Cerebrospinal fluid pressure (~10 mmHg)

# Lamina cribrosa compliance derivation (Winkler elastic foundation)
E_LC_MPa = 0.3  # Average Young's modulus of LC (~0.3 MPa for healthy, higher for glaucomatous)
t_LC_mm = 0.3   # Average thickness of LC (~0.3 mm / 300 micrometers)
k_LC = Constant(E_LC_MPa / t_LC_mm) # Calculated spring constant (MPa/mm)

zero = Constant((0.0, 0.0))

axis_id = ID_SYMMETRY_AXIS
bc_sym_axis = DirichletBC(V.sub(0), Constant(0.0), boundaries, axis_id)

bottom_robin_id = ID_LC_FLOOR
# This is now a Robin BC (compliant spring), replacing the old rigid roller BC. See Pi_total for formulation.

u_outer_rad = Constant(0.0) # Scleral expansion tracker

# Apply radial expansion and prevent vertical movement on the Sclera Wall
bc_sclera_wall_x = DirichletBC(V.sub(0), u_outer_rad, boundaries, ID_SCLERA_WALL)
bc_sclera_wall_y = DirichletBC(V.sub(1), Constant(0.0), boundaries, ID_SCLERA_WALL)

# Prevent vertical movement on the Retina Boundary Roller (freely moves horizontally)
bc_retina_roller_y = DirichletBC(V.sub(1), Constant(0.0), boundaries, ID_RETINA_BOUNDARY)

# Collect all boundary conditions
all_bcs = [bc_sym_axis, bc_sclera_wall_x, bc_sclera_wall_y, bc_retina_roller_y]

nu = 0.499                                  # Poisson's ratio (near incompressible)
mu_val = 56.0 / 1000                        # Float value for Python calculations, convert kPa to MPa
E_val = 2 * mu_val * (1 + nu)               # Back cal the Young's mod
lmbda_val = (E_val * nu) / ((1 + nu) * (1 - 2 * nu)) # Calculate lmbda as a float

V_DG0_scalar = FunctionSpace(mesh, "DG", 0)
V_DG0_tensor = TensorFunctionSpace(mesh, "DG", 0, shape=(3, 3))

mu = Function(V_DG0_scalar)
mu.vector()[:] = mu_val
lmbda = Function(V_DG0_scalar)
lmbda.vector()[:] = lmbda_val

#healthy baseline expression of GFAP and ECM
rho_gfap = Function(V_DG0_scalar)
rho_gfap.vector()[:] = 1.0
rho_ecm = Function(V_DG0_scalar)
rho_ecm.vector()[:] = 1.0

u = Function(V)
u.rename("Displacement", "Displacement")
du = TrialFunction(V)
v = TestFunction(V)

I3 = Identity(3)

# Helper function to construct 3D axisymmetric deformation gradient
def get_F3D(u):
    r = SpatialCoordinate(mesh)[0]
    r_safe = r + 1e-14 # small epsilon to prevent division by zero at symmetry axis
    u_r, u_z = u[0], u[1]
    grad_u = grad(u)
    return as_tensor([
        [1 + grad_u[0, 0], 0, grad_u[0, 1]],
        [0, 1 + u_r / r_safe, 0],
        [grad_u[1, 0], 0, 1 + grad_u[1, 1]]
    ])

vm_func = Function(V_DG0_scalar)
vm_func.rename("von_Mises_Stress", "von_Mises_Stress")
E3_func = Function(V_DG0_scalar)
E3_func.rename("Min_Principal_Strain", "Min_Principal_Strain")
hoop_func = Function(V_DG0_scalar)
hoop_func.rename("Hoop_Stress", "Hoop_Stress")
E1_func = Function(V_DG0_scalar)
E1_func.rename("Max_Principal_Strain", "Max_Principal_Strain")

# Define a custom measure with fixed quadrature degree for consistency
dx_c = Measure("dx", domain=mesh, metadata={"quadrature_degree": 4})
r_coord = SpatialCoordinate(mesh)[0]
r_fac = 2 * np.pi * r_coord

# ======================================================
# 4. Loading Phase (step load)
# ======================================================
load_steps = 100
pressure_id = ID_CUP_SURFACE_IOP # ID for IOP boundary
load_factor = Constant(0.0)
ds_measure = Measure("ds", domain=mesh, subdomain_data=boundaries)
n = FacetNormal(mesh)

F = get_F3D(u)
C = F.T * F
Ic = tr(C)
J = det(F)
psi = (mu / 2) * (Ic - 3) - mu * ln(J) + (lmbda / 2) * (ln(J))**2   #neo hookean equation

# Potential energy with normal pressure and bottom Robin BC
Pi_total = (psi * r_fac * dx_c 
            + load_factor * iop_target_MPa * dot(n, u) * r_fac * ds_measure(pressure_id)
            + load_factor * csfp_target_MPa * dot(n, u) * r_fac * ds_measure(bottom_robin_id)
            + (k_LC / 2.0) * (u[1]**2) * r_fac * ds_measure(bottom_robin_id))
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

# UFL for exact Cauchy stress during loading (Fg = I)
F_var_L = variable(get_F3D(u))
C_L = F_var_L.T * F_var_L
Ic_L = tr(C_L)
J_L = det(F_var_L)
Psi_L = (mu / 2) * (Ic_L - 3) - mu * ln(J_L) + (lmbda / 2) * (ln(J_L))**2
P_L = diff(Psi_L, F_var_L)
sigma_L = (1.0 / J_L) * P_L * F_var_L.T
sigma_dev_L = sigma_L - (1.0 / 3.0) * tr(sigma_L) * I3
vm_expr_L = sqrt((3.0 / 2.0) * inner(sigma_dev_L, sigma_dev_L))

# Calculate Minimum Principal Stress (Max Compressive)
sig_center_L = (sigma_L[0, 0] + sigma_L[2, 2]) / 2.0
sig_radius_L = sqrt(((sigma_L[0, 0] - sigma_L[2, 2]) / 2.0)**2 + sigma_L[0, 2]**2)
sig_3_L = conditional(lt(sig_center_L - sig_radius_L, sigma_L[1, 1]), sig_center_L - sig_radius_L, sigma_L[1, 1])
p_expr_L = sig_3_L
hoop_expr_L = sigma_L[1, 1]

# Calculate Maximum and Minimum Principal Green-Lagrange Strains
E_L = 0.5 * (C_L - I3)
E_center_L = (E_L[0, 0] + E_L[2, 2]) / 2.0
E_radius_L = sqrt(((E_L[0, 0] - E_L[2, 2]) / 2.0)**2 + E_L[0, 2]**2)
E_1_L = conditional(gt(E_center_L + E_radius_L, E_L[1, 1]), E_center_L + E_radius_L, E_L[1, 1])
E_3_L = conditional(lt(E_center_L - E_radius_L, E_L[1, 1]), E_center_L - E_radius_L, E_L[1, 1])

# ======================================================
# Pre-calculate Maximum Dynamic Scleral Expansion
# ======================================================
R_eye = 12.0 # mm
t_sclera = 1.0 # mm
E_sclera = 3.0 # MPa

# 1. Nominal Scleral Tension (Laplace's Law) at maximum target pressure
sigma_0_max = (iop_target_MPa * R_eye) / (2.0 * t_sclera)
# 2. Stress Concentration (Kirsch's Equation) at the edge of the hole i can use 2 directly
sigma_theta_max = 2.0 * sigma_0_max
# 3. Circumferential Stretch (Hooke's Law)
epsilon_theta_max = (1.0 / E_sclera) * sigma_theta_max
# 4. Total Maximum Expansion
canal_radius = 0.8 # Inner radius of the scleral canal
Delta_r_max = canal_radius * epsilon_theta_max

print("=== Loading Phase ===")
for step in range(1, load_steps + 1):
    factor = step / load_steps
    load_factor.assign(factor) # load_factor now scales from 0 to 1
    
    u_outer_rad.assign(Delta_r_max * factor)
    solve(F_res == 0, u, bcs=all_bcs, J=J_res, solver_parameters=solver_params)
    
    vm_func.assign(project(vm_expr_L, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    E3_func.assign(project(E_3_L, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    hoop_func.assign(project(hoop_expr_L, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    E1_func.assign(project(E_1_L, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    
    xdmf_load.write(u, float(step))
    xdmf_load.write(vm_func, float(step))
    xdmf_load.write(E3_func, float(step))
    xdmf_load.write(hoop_func, float(step))
    xdmf_load.write(E1_func, float(step))
xdmf_load.close()
print("✅ Loading phase finished and saved.")

# Capture the deposition state (Reference Configuration for Scar Tissue)
F_deposit = project(get_F3D(u), V_DG0_tensor, form_compiler_parameters={"quadrature_degree": 2})

mu_scar_func = Function(V_DG0_scalar)
mu_scar_func.vector()[:] = 0.0
lmbda_scar_func = Function(V_DG0_scalar)
lmbda_scar_func.vector()[:] = 0.0

# ======================================================
# 5. Remodeling Phase (maintain load)
# ======================================================
growth_steps = 500
dt_growth = 0.1
growth_rate = 0.05

theta_func = Function(V_DG0_scalar)
theta_func.vector()[:] = 1.0

# --- Solve Laplace PDE for Axon Orientation ---
# We use a potential field to naturally guide axons from the retina (source) to the LC (sink)
V_phi = FunctionSpace(mesh, "CG", 2)
phi = TrialFunction(V_phi)
v_phi = TestFunction(V_phi)

# Axisymmetric weak form of Laplace equation: div(grad(phi)) = 0
a_phi = inner(grad(phi), grad(v_phi)) * r_fac * dx_c
L_phi = Constant(0.0) * v_phi * r_fac * dx_c

# Boundary conditions: High potential at Retina, Low potential at LC
# Natural boundary conditions (zero flux) on other walls force flow perfectly parallel to them
bc_source = DirichletBC(V_phi, Constant(1.0), boundaries, ID_RETINA_BOUNDARY)
bc_sink = DirichletBC(V_phi, Constant(0.0), boundaries, ID_LC_FLOOR)

phi_sol = Function(V_phi)
solve(a_phi == L_phi, phi_sol, [bc_source, bc_sink])

# The axons flow down the gradient (from 1.0 to 0.0)
grad_phi = grad(phi_sol)
grad_norm = sqrt(inner(grad_phi, grad_phi)) + 1e-14 # Epsilon prevents division by zero
n0_2d = -grad_phi / grad_norm

# Map to 3D axisymmetric coordinates (r, theta, z) -> (n_r, 0, n_z)
n0 = as_vector((n0_2d[0], 0.0, n0_2d[1]))

# --- Export Axon Orientation for ParaView Visualization ---
V_DG0_vec = VectorFunctionSpace(mesh, "DG", 0)
n0_func = project(n0_2d, V_DG0_vec, form_compiler_parameters={"quadrature_degree": 2})
n0_func.rename("Axon_Orientation", "Reference Axon Direction")

with XDMFFile(mesh.mpi_comm(), "axon_orientation.xdmf") as xdmf_n0:
    xdmf_n0.write(n0_func)

# Keep load at final value
load_factor.assign(1.0) # Keep load at full target_load (1.0 * iop_target_MPa)

u_outer_rad.assign(Delta_r_max)

xdmf_growth = XDMFFile(mesh.mpi_comm(), "u_growth_phase.xdmf")
xdmf_growth.parameters["flush_output"] = True
xdmf_growth.parameters["functions_share_mesh"] = True

xdmf_theta = XDMFFile(mesh.mpi_comm(), "theta_growth_phase.xdmf")
xdmf_theta.parameters["flush_output"] = True
xdmf_theta.parameters["functions_share_mesh"] = True

# --- Define variational problem once outside the loop ---
F_current_sym = get_F3D(u)
Fg_ufl_sym = sqrt(theta_func) * I3 + (1.0 - sqrt(theta_func)) * outer(n0, n0)

# 1. Neural Component
Fe_neural_sym = F_current_sym * inv(Fg_ufl_sym)
C_neural_sym = Fe_neural_sym.T * Fe_neural_sym
Ic_neural_sym = tr(C_neural_sym)
J_neural_sym = det(Fe_neural_sym)
psi_neural = (mu_val/2)*(Ic_neural_sym - 3) - mu_val*ln(J_neural_sym) + (lmbda_val/2)*(ln(J_neural_sym))**2 

# 2. Scar Component (Deposited at F_deposit, generates 0 stress unless deformed further)
Fe_scar_sym = F_current_sym * inv(F_deposit)
C_scar_sym = Fe_scar_sym.T * Fe_scar_sym
Ic_scar_sym = tr(C_scar_sym)
J_scar_sym = det(Fe_scar_sym)
psi_scar = (mu_scar_func/2)*(Ic_scar_sym - 3) - mu_scar_func*ln(J_scar_sym) + (lmbda_scar_func/2)*(ln(J_scar_sym))**2 

psi_growth = psi_neural + psi_scar
Pi_total_growth = (psi_growth * r_fac * dx_c 
                   + load_factor * iop_target_MPa * dot(n, u) * r_fac * ds_measure(pressure_id)
                   + load_factor * csfp_target_MPa * dot(n, u) * r_fac * ds_measure(bottom_robin_id)
                   + (k_LC / 2.0) * (u[1]**2) * r_fac * ds_measure(bottom_robin_id))
F_res_growth = derivative(Pi_total_growth, u, v)
J_res_growth = derivative(F_res_growth, u, du)

print("=== Remodeling Phase (Atrophy/Shrinkage + Stiffening) ===")

# Open a log file to record remodeling status ("w" mode ensures it overwrites any previous file)
log_file = open("remodeling_log.txt", "w")
log_file.write("Step,Atrophy Rate,Avg mu,Max mu,Std mu\n")

# Biological ODE Parameters
strain_homeo = 1.8e-2  # Dimensionless (Homeostatic Green-Lagrange strain baseline: 0.28% stretch)
k_gfap = 0.0010      # GFAP synthesis rate (calibrated for 63 kPa max)
d_gfap = 0.20        # GFAP degradation rate (calibrated for step 250 plateau)
k_ecm = 0.0010       # ECM synthesis rate (calibrated for 63 kPa max)
d_ecm = 0.20         # ECM degradation rate (calibrated for step 250 plateau)
c1 = 0.0030          # Stiffness contribution of GFAP
c2 = 0.0070          # Stiffness contribution of ECM

# Atrophy/Shrinkage Parameters
theta_min = 0.3      # Maximum allowed physiological shrinkage (30%)
gamma_atrophy = 2.0  # Exponent for smooth saturation (1=linear, 2=quadratic)

atrophy_field = Function(V_DG0_scalar)
mu_scar_func_old = Function(V_DG0_scalar)
delta_mu_scar_func = Function(V_DG0_scalar)

for gstep in range(1, growth_steps + 1):
    # --- 1. Element-wise Stress Extraction ---
    F_var = variable(get_F3D(u))
    Fg_ufl_var = sqrt(theta_func) * I3 + (1.0 - sqrt(theta_func)) * outer(n0, n0)
    
    # Neural Stress
    Fe_neural_var = F_var * inv(Fg_ufl_var)
    C_neural_var = Fe_neural_var.T * Fe_neural_var
    Ic_neural_var = tr(C_neural_var)
    J_neural_var = det(Fe_neural_var)
    psi_neural_var = (mu_val / 2) * (Ic_neural_var - 3) - mu_val * ln(J_neural_var) + (lmbda_val / 2) * (ln(J_neural_var))**2
    
    # Scar Stress
    Fe_scar_var = F_var * inv(F_deposit)
    C_scar_var = Fe_scar_var.T * Fe_scar_var
    Ic_scar_var = tr(C_scar_var)
    J_scar_var = det(Fe_scar_var)
    psi_scar_var = (mu_scar_func / 2) * (Ic_scar_var - 3) - mu_scar_func * ln(J_scar_var) + (lmbda_scar_func / 2) * (ln(J_scar_var))**2
    
    # Total Combined Stress
    Psi_tot_var = psi_neural_var + psi_scar_var
    P_exact = diff(Psi_tot_var, F_var)
    J_tot = det(F_var)
    sigma_exact = (1.0 / J_tot) * P_exact * F_var.T
    
    # Calculate Principal Strains (E1 for stiffening, E3 for capillary collapse)
    E_neural_var = 0.5 * (C_neural_var - I3)
    E_center = (E_neural_var[0, 0] + E_neural_var[2, 2]) / 2.0
    E_radius = sqrt(((E_neural_var[0, 0] - E_neural_var[2, 2]) / 2.0)**2 + E_neural_var[0, 2]**2)
    E_1 = conditional(gt(E_center + E_radius, E_neural_var[1, 1]), E_center + E_radius, E_neural_var[1, 1])
    E_3 = conditional(lt(E_center - E_radius, E_neural_var[1, 1]), E_center - E_radius, E_neural_var[1, 1])
    E3_func.assign(project(E_3, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    
    # Exact von Mises Stress
    sigma_dev_exact = sigma_exact - (1.0 / 3.0) * tr(sigma_exact) * I3
    vm_expr_G = sqrt((3.0 / 2.0) * inner(sigma_dev_exact, sigma_dev_exact))
    vm_func.assign(project(vm_expr_G, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    
    # Project Hoop Stress
    hoop_func.assign(project(sigma_exact[1, 1], V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    
    # Project Max Principal Strain for ParaView
    E1_func.assign(project(E_1, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    
    # Project Maximum Principal Strain for Stiffening Trigger
    strain_mag_expr = E_1
    strain_mag_func = project(strain_mag_expr, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2})
    
    # Project Growth Tensor Volume (Death Trigger)
    J_g_expr = det(Fg_ufl_var)
    J_g_func = project(J_g_expr, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2})
    
    # Project Elastic Stretch (Mechanical Rupture Death Trigger)
    # Ic_e measures the pure elastic distortion/stretch of the cell
    Ic_e_func = project(Ic_neural_var, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2})
    
    # Save previous scar stiffness for F_deposit weighting
    mu_scar_func_old.assign(mu_scar_func)
    
    # --- 2. Fast NumPy Biological Math (Atrophy + Stiffening) ---
    E3_array = E3_func.vector().get_local()
    strain_mag_array = strain_mag_func.vector().get_local()
    mu_array = mu.vector().get_local()
    lmbda_array = lmbda.vector().get_local()
    rho_gfap_arr = rho_gfap.vector().get_local()
    rho_ecm_arr = rho_ecm.vector().get_local()
    J_g_array = J_g_func.vector().get_local()
    Ic_e_array = Ic_e_func.vector().get_local()
    theta_array = theta_func.vector().get_local()
    mu_scar_arr_old = mu_scar_func_old.vector().get_local()
    
    # Mechanobiological Stiffening (ODE integration)
    Phi_arr = np.maximum(0.0, (strain_mag_array - strain_homeo) / strain_homeo)
    rho_gfap_arr += dt_growth * (k_gfap * Phi_arr - d_gfap * (rho_gfap_arr - 1.0))
    rho_ecm_arr += dt_growth * (k_ecm * Phi_arr - d_ecm * (rho_ecm_arr - 1.0))
    
    # Calculate Scar Stiffness Components
    mu_scar_arr_new = c1 * (rho_gfap_arr - 1.0) + c2 * (rho_ecm_arr - 1.0)
    delta_mu_scar_arr = mu_scar_arr_new - mu_scar_arr_old
    lmbda_scar_arr = lmbda_val * (mu_scar_arr_new / mu_val)
    
    # Total Apparent Stiffness (for logging/ParaView visualization only)
    mu_array = mu_val + mu_scar_arr_new
    
    # Calculate local blood flow reduction based on Capillary Squishing (Poiseuille's Law)
    # Green-Lagrange Strain to Stretch Ratio squared: lambda_3^2 = 2*E3 + 1
    # Poiseuille Flow Q is proportional to r^4 -> Q_local = (lambda_3^2)^2
    # We clip to prevent negative stretch values in extreme non-physical compressions
    lambda3_sq = np.clip(2.0 * E3_array + 1.0, 0.0, None)
    Q_local_array = lambda3_sq**2
    
    # Compute local atrophy trigger 
    # E3_crit = -0.032 -> lambda3_sq = 0.936 -> Q_crit = 0.876096
    atrophy_array = np.where(Q_local_array < 0.876096, 0.876096 - Q_local_array, 0.0)
    
    # --- 3. Update FEniCS Fields ---
    mu.vector().set_local(mu_array)
    mu.vector().apply("insert")
    rho_gfap.vector().set_local(rho_gfap_arr)
    rho_gfap.vector().apply("insert")
    rho_ecm.vector().set_local(rho_ecm_arr)
    rho_ecm.vector().apply("insert")
    mu_scar_func.vector().set_local(mu_scar_arr_new)
    mu_scar_func.vector().apply("insert")
    delta_mu_scar_func.vector().set_local(delta_mu_scar_arr)
    delta_mu_scar_func.vector().apply("insert")
    lmbda_scar_func.vector().set_local(lmbda_scar_arr)
    lmbda_scar_func.vector().apply("insert")

    atrophy_field.vector().set_local(atrophy_array)
    atrophy_field.vector().apply("insert")
    
    # --- Atrophy/Shrinkage increment ---
    # Apply Kuhl's smooth saturation function (adapted for atrophy)
    # k_theta smoothly slows down the atrophy rate as theta approaches theta_min
    k_theta = growth_rate * ((theta_array - theta_min) / (1.0 - theta_min))**gamma_atrophy
    k_theta = np.maximum(k_theta, 0.0) # Prevent negative rates
    
    theta_array -= k_theta * atrophy_array * dt_growth * theta_array
    theta_array = np.clip(theta_array, theta_min, 1.0) # Retained purely as a numerical safety net
    theta_func.vector().set_local(theta_array)
    theta_func.vector().apply("insert")

    # --- Constrained Mixture: Update Reference Configuration for New Scar Tissue ---
    # New ECM is deposited at the current stretched state. The new effective reference configuration 
    # is a stiffness-weighted average of the old reference and the current deformation.
    mu_scar_safe = conditional(lt(mu_scar_func, 1e-12), 1e-12, mu_scar_func)
    F_deposit_expr = conditional(gt(delta_mu_scar_func, 0.0),
                                 (mu_scar_func_old * F_deposit + delta_mu_scar_func * get_F3D(u)) / mu_scar_safe,
                                 F_deposit)
    F_deposit.assign(project(F_deposit_expr, V_DG0_tensor, form_compiler_parameters={"quadrature_degree": 2}))

    # --- Solve ---
    # The solver 're-reads' the updated theta values inside the symbolic F_res_growth form
    solve(F_res_growth == 0, u, bcs=all_bcs, J=J_res_growth, solver_parameters=solver_params)

    # --- Save ---
    if gstep % 10 == 0:
        avg_atrophy = np.mean(atrophy_array)
        avg_mu = np.mean(mu_array)
        max_mu = np.max(mu_array)
        std_mu = np.std(mu_array)
        print(f"[Remodeling Step {gstep}] Atrophy Rate = {avg_atrophy:.4e}")
        print(f"                     Avg mu = {avg_mu:.4f} (Max: {max_mu:.4f}, StdDev: {std_mu:.4f})")
        print(f"                     Max I_c (Distortion) = {np.max(Ic_e_array):.4f}")
        
        # Write status to log file
        log_line = f"{gstep},{avg_atrophy:.4e},{avg_mu:.4f},{max_mu:.4f},{std_mu:.4f}\n"
        log_file.write(log_line)

    xdmf_growth.write(u, float(gstep))
    xdmf_growth.write(vm_func, float(gstep))
    xdmf_growth.write(E3_func, float(gstep))
    xdmf_growth.write(hoop_func, float(gstep))
    xdmf_growth.write(E1_func, float(gstep))
    if gstep % 20 == 0:
        xdmf_theta.write(theta_func, float(gstep))

xdmf_growth.close()
xdmf_theta.close()
log_file.close()
print("✅ Remodeling phase finished and saved.")

# ======================================================
# 6. Unloading Phase (gradual release from remodeled state)
# ======================================================
unload_steps = 100
xdmf_unload = XDMFFile(mesh.mpi_comm(), "u_unloading_phase.xdmf")
xdmf_unload.parameters["flush_output"] = True
xdmf_unload.parameters["functions_share_mesh"] = True

# --- Define variational problem once for unloading ---
F_current_un = get_F3D(u)
Fg_ufl_un = sqrt(theta_func) * I3 + (1.0 - sqrt(theta_func)) * outer(n0, n0)

Fe_neural_un = F_current_un * inv(Fg_ufl_un)
C_neural_un = Fe_neural_un.T * Fe_neural_un
Ic_neural_un = tr(C_neural_un)
J_neural_un = det(Fe_neural_un)
psi_neural_un = (mu_val / 2) * (Ic_neural_un - 3) - mu_val * ln(J_neural_un) + (lmbda_val / 2) * (ln(J_neural_un))**2

Fe_scar_un = F_current_un * inv(F_deposit)
C_scar_un = Fe_scar_un.T * Fe_scar_un
Ic_scar_un = tr(C_scar_un)
J_scar_un = det(Fe_scar_un)
psi_scar_un = (mu_scar_func / 2) * (Ic_scar_un - 3) - mu_scar_func * ln(J_scar_un) + (lmbda_scar_func / 2) * (ln(J_scar_un))**2

psi_un = psi_neural_un + psi_scar_un
Pi_total_unload = (psi_un * r_fac * dx_c 
                   + load_factor * iop_target_MPa * dot(n, u) * r_fac * ds_measure(pressure_id)
                   + load_factor * csfp_target_MPa * dot(n, u) * r_fac * ds_measure(bottom_robin_id)
                   + (k_LC / 2.0) * (u[1]**2) * r_fac * ds_measure(bottom_robin_id))
F_res_unload = derivative(Pi_total_unload, u, v)
J_res_unload = derivative(F_res_unload, u, du)

# UFL for exact Cauchy stress during unloading
F_var_U = variable(get_F3D(u))
Fg_ufl_U = sqrt(theta_func) * I3 + (1.0 - sqrt(theta_func)) * outer(n0, n0)

Fe_neural_U = F_var_U * inv(Fg_ufl_U)
C_neural_U = Fe_neural_U.T * Fe_neural_U
Ic_neural_U = tr(C_neural_U)
J_neural_U = det(Fe_neural_U)
psi_neural_U = (mu_val / 2) * (Ic_neural_U - 3) - mu_val * ln(J_neural_U) + (lmbda_val / 2) * (ln(J_neural_U))**2

Fe_scar_U = F_var_U * inv(F_deposit)
C_scar_U = Fe_scar_U.T * Fe_scar_U
Ic_scar_U = tr(C_scar_U)
J_scar_U = det(Fe_scar_U)
psi_scar_U = (mu_scar_func / 2) * (Ic_scar_U - 3) - mu_scar_func * ln(J_scar_U) + (lmbda_scar_func / 2) * (ln(J_scar_U))**2

Psi_U = psi_neural_U + psi_scar_U
P_exact_U = diff(Psi_U, F_var_U)
J_tot_U = det(F_var_U)
sigma_exact_U = (1.0 / J_tot_U) * P_exact_U * F_var_U.T
sigma_dev_U = sigma_exact_U - (1.0 / 3.0) * tr(sigma_exact_U) * I3
vm_expr_U = sqrt((3.0 / 2.0) * inner(sigma_dev_U, sigma_dev_U))

hoop_expr_U = sigma_exact_U[1, 1]

# Calculate Maximum and Minimum Principal Green-Lagrange Strains
E_neural_U = 0.5 * (C_neural_U - I3)
E_center_U = (E_neural_U[0, 0] + E_neural_U[2, 2]) / 2.0
E_radius_U = sqrt(((E_neural_U[0, 0] - E_neural_U[2, 2]) / 2.0)**2 + E_neural_U[0, 2]**2)
E_1_U = conditional(gt(E_center_U + E_radius_U, E_neural_U[1, 1]), E_center_U + E_radius_U, E_neural_U[1, 1])
E_3_U = conditional(lt(E_center_U - E_radius_U, E_neural_U[1, 1]), E_center_U - E_radius_U, E_neural_U[1, 1])

print("=== Unloading Phase ===")
for step in range(1, unload_steps + 1):
    # Gradually reduce load and scleral expansion
    factor = 1.0 * (1 - step / unload_steps)
    load_factor.assign(factor) # Scale from 1.0 down to 0.0
    
    u_outer_rad.assign(Delta_r_max * factor)

    # --- Solve equilibrium ---
    solve(F_res_unload == 0, u, bcs=all_bcs, J=J_res_unload, solver_parameters=solver_params)
    
    vm_func.assign(project(vm_expr_U, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    E3_func.assign(project(E_3_U, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    hoop_func.assign(project(hoop_expr_U, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    E1_func.assign(project(E_1_U, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    
    xdmf_unload.write(u, float(step))
    xdmf_unload.write(vm_func, float(step))
    xdmf_unload.write(E3_func, float(step))
    xdmf_unload.write(hoop_func, float(step))
    xdmf_unload.write(E1_func, float(step))

xdmf_unload.close()
print("✅ Unloading phase finished and saved.")
