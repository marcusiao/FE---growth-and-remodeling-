from dolfin import *
import numpy as np
import os

# ======================================================
# 1. Mesh Generation (2D Prelaminar Tissue)
# ======================================================
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Define domain resolution and extents
nx, ny = 60, 40
x_min, x_max = 0.0, 2.0
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

class SymmetryAxis(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], x_min)

class Top(SubDomain):
    def inside(self, x, on_boundary):
        # Mark the top curve by excluding the vertical side walls (at x_min and x_max)
        return on_boundary and x[1] > 0.1 and not (near(x[0], x_min) or near(x[0], x_max))

class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], y_min)

class OuterSide(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], x_max)

SymmetryAxis().mark(boundaries, 17) # axis_id
Top().mark(boundaries, 14)    # pressure_id
Bottom().mark(boundaries, 16) # bottom_roller_id
OuterSide().mark(boundaries, 18) # outer_side_id

V = VectorFunctionSpace(mesh, "CG", 2)  #2 mean quadratic

# ======================================================
# 3. Boundary conditions and material parameters
# ======================================================
# Consistent Unit System (Standard for Biomechanics):
# Length: mm | Force: N | Stress/Pressure: MPa (N/mm^2)
# 1 mmHg = 0.0001333 MPa

MPa_TO_mmHg = 1.0 / 0.0001333 # Conversion factor

# Target IOP (will be scaled by load_factor)
iop_target_MPa = 40 * 0.0001333  # Increased from 40 to 150 mmHg to increase tissue stress
lc_target_MPa = 10 * 0.0001333   # LC pressure (not being used)


zero = Constant((0.0, 0.0))

axis_id = 17
bc_sym_axis = DirichletBC(V.sub(0), Constant(0.0), boundaries, axis_id)

bottom_roller_id = 16
bc_bottom_roller = DirichletBC(V.sub(1), Constant(0.0), boundaries, bottom_roller_id)

outer_side_id = 18
u_outer_rad = Constant(0.0) # Scleral expansion tracker
bc_outer_roller = DirichletBC(V.sub(0), u_outer_rad, boundaries, outer_side_id)

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

d = mesh.geometry().dim()
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
p_func = Function(V_DG0_scalar)
p_func.rename("Min_Principal_Stress", "Min_Principal_Stress")
hoop_func = Function(V_DG0_scalar)
hoop_func.rename("Hoop_Stress", "Hoop_Stress")
sig1_func = Function(V_DG0_scalar)
sig1_func.rename("Max_Principal_Stress", "Max_Principal_Stress")

# Define a custom measure with fixed quadrature degree for consistency
dx_c = Measure("dx", domain=mesh, metadata={"quadrature_degree": 4})
r_coord = SpatialCoordinate(mesh)[0]
r_fac = 2 * np.pi * r_coord

# ======================================================
# 4. Loading Phase (step load)
# ======================================================
load_steps = 100
pressure_id = 14 # ID for IOP boundary
load_factor = Constant(0.0)
ds_measure = Measure("ds", domain=mesh, subdomain_data=boundaries)
n = FacetNormal(mesh)

F = get_F3D(u)
C = F.T * F
Ic = tr(C)
J = det(F)
psi = (mu / 2) * (Ic - 3) - mu * ln(J) + (lmbda / 2) * (ln(J))**2   #neo hookean equation

# Potential energy with normal pressure: Pi = Psi - (-P * n . u) = Psi + P * n . u
Pi_total = psi * r_fac * dx_c + load_factor * iop_target_MPa * dot(n, u) * r_fac * ds_measure(pressure_id)
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
sig_1_L = conditional(gt(sig_center_L + sig_radius_L, sigma_L[1, 1]), sig_center_L + sig_radius_L, sigma_L[1, 1])

print("=== Loading Phase ===")
for step in range(1, load_steps + 1):
    factor = step / load_steps
    load_factor.assign(factor) # load_factor now scales from 0 to 1
    u_outer_rad.assign(0.05 * x_max * factor) # 5% radial expansion ramp
    solve(F_res == 0, u, bcs=[bc_sym_axis, bc_bottom_roller, bc_outer_roller], J=J_res, solver_parameters=solver_params)
    
    vm_func.assign(project(vm_expr_L, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    p_func.assign(project(p_expr_L, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    hoop_func.assign(project(hoop_expr_L, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    sig1_func.assign(project(sig_1_L, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    
    xdmf_load.write(u, float(step))
    xdmf_load.write(vm_func, float(step))
    xdmf_load.write(p_func, float(step))
    xdmf_load.write(hoop_func, float(step))
    xdmf_load.write(sig1_func, float(step))
xdmf_load.close()
print("✅ Loading phase finished and saved.")

# ======================================================
# 5. Remodeling Phase (maintain load)
# ======================================================
growth_steps = 500
dt_growth = 0.1
growth_rate = 0.05

theta_func = Function(V_DG0_scalar)
theta_func.vector()[:] = 1.0

# Define spatially varying axon orientation (horizontal at sides, turning vertical at the center canal)
x_coords = SpatialCoordinate(mesh)
# In axisymmetry, r = x_coords[0] is always >= 0.
angle_from_vertical = conditional(lt(x_coords[0], 0.5), np.pi * x_coords[0], np.pi/2)
# Create the 2D vector by rotating a vertical vector (0,1) by this angle.
n0_2d = as_vector((sin(angle_from_vertical), cos(angle_from_vertical)))
# Map to 3D axisymmetric coordinates (r, theta, z) -> (n_r, 0, n_z)
n0 = as_vector((n0_2d[0], 0.0, n0_2d[1]))

# Keep load at final value
load_factor.assign(1.0) # Keep load at full target_load (1.0 * iop_target_MPa)
u_outer_rad.assign(0.05 * x_max) # Keep scleral expansion at 5%

xdmf_growth = XDMFFile(mesh.mpi_comm(), "u_growth_phase.xdmf")
xdmf_growth.parameters["flush_output"] = True
xdmf_growth.parameters["functions_share_mesh"] = True

xdmf_theta = XDMFFile(mesh.mpi_comm(), "theta_growth_phase.xdmf")
xdmf_theta.parameters["flush_output"] = True
xdmf_theta.parameters["functions_share_mesh"] = True

# --- Define variational problem once outside the loop ---
F_current_sym = get_F3D(u)
Fg_ufl_sym = sqrt(theta_func) * I3 + (1.0 - sqrt(theta_func)) * outer(n0, n0)
Fe_sym = F_current_sym * inv(Fg_ufl_sym)  # This is the elastic part Fe from F = Fe * Fg
C_sym = Fe_sym.T * Fe_sym
Ic_sym = tr(C_sym)
J_sym = det(Fe_sym)
psi_growth = (mu/2)*(Ic_sym - 3) - mu*ln(J_sym) + (lmbda/2)*(ln(J_sym))**2 # Use current mu, lmbda
Pi_total_growth = psi_growth * r_fac * dx_c + load_factor * iop_target_MPa * dot(n, u) * r_fac * ds_measure(pressure_id)
F_res_growth = derivative(Pi_total_growth, u, v)
J_res_growth = derivative(F_res_growth, u, du)

print("=== Remodeling Phase (Atrophy/Shrinkage + Stiffening) ===")
current_iop_mmHg = iop_target_MPa * MPa_TO_mmHg

# Open a log file to record remodeling status ("w" mode ensures it overwrites any previous file)
log_file = open("remodeling_log.txt", "w")
log_file.write("Step,Atrophy Rate,Avg mu,Max mu,Std mu\n")

# Biological ODE Parameters
sigma_homeo = 0.00041  # MPa (Homeostatic von Mises stress baseline)
k_gfap = 0.005       # GFAP synthesis rate
d_gfap = 0.05        # GFAP degradation rate (raised to create a natural biological ceiling)
k_ecm = 0.004        # ECM synthesis rate (calibrated for a natural plateau)
d_ecm = 0.04         # ECM degradation rate (calibrated for a natural plateau)
c1 = 0.0042          # Stiffness contribution of GFAP (10 kPa per unit density)
c2 = 0.0084         # Stiffness contribution of ECM (100 kPa per unit density - Dense Collagen)

atrophy_field = Function(V_DG0_scalar)

for gstep in range(1, growth_steps + 1):
    # --- 1. Element-wise Stress Extraction ---
    F_var = variable(get_F3D(u))
    Fg_ufl_var = sqrt(theta_func) * I3 + (1.0 - sqrt(theta_func)) * outer(n0, n0)
    Fe_var = F_var * inv(Fg_ufl_var)
    
    # Compute exact Cauchy stress analytically using UFL
    C_e = Fe_var.T * Fe_var
    Ic_e = tr(C_e)
    J_e = det(Fe_var)
    J_tot = det(F_var)
    
    Psi = (mu / 2) * (Ic_e - 3) - mu * ln(J_e) + (lmbda / 2) * (ln(J_e))**2
    P_exact = diff(Psi, F_var)
    sigma_exact = (1.0 / J_tot) * P_exact * F_var.T
    
    # Calculate Minimum Principal Stress (Max Compressive) to replace Hydrostatic Pressure
    sig_center = (sigma_exact[0, 0] + sigma_exact[2, 2]) / 2.0
    sig_radius = sqrt(((sigma_exact[0, 0] - sigma_exact[2, 2]) / 2.0)**2 + sigma_exact[0, 2]**2)
    sig_3 = conditional(lt(sig_center - sig_radius, sigma_exact[1, 1]), sig_center - sig_radius, sigma_exact[1, 1])
    p_tissue_expr = sig_3
    
    # Calculate Maximum Principal Stress (Max Tensile) for Astrocyte Stiffening Trigger
    sig_1 = conditional(gt(sig_center + sig_radius, sigma_exact[1, 1]), sig_center + sig_radius, sigma_exact[1, 1])
    
    p_func.assign(project(p_tissue_expr, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    
    # Exact von Mises Stress
    sigma_dev_exact = sigma_exact - (1.0 / 3.0) * tr(sigma_exact) * I3
    vm_expr_G = sqrt((3.0 / 2.0) * inner(sigma_dev_exact, sigma_dev_exact))
    vm_func.assign(project(vm_expr_G, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    
    # Project Hoop Stress
    hoop_func.assign(project(sigma_exact[1, 1], V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    
    # Project Max Principal Stress for ParaView
    sig1_func.assign(project(sig_1, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    
    # Project Maximum Principal Stress for Stiffening Trigger (Replaces von Mises)
    sigma_mag_expr = sig_1
    sigma_mag_func = project(sigma_mag_expr, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2})
    
    # Project Growth Tensor Volume (Death Trigger)
    J_g_expr = det(Fg_ufl_var)
    J_g_func = project(J_g_expr, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2})
    
    # Project Elastic Stretch (Mechanical Rupture Death Trigger)
    # Ic_e measures the pure elastic distortion/stretch of the cell
    Ic_e_func = project(Ic_e, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2})
    
    # --- 2. Fast NumPy Biological Math (Atrophy + Stiffening) ---
    p_tissue_array = p_func.vector().get_local()
    sigma_mag_array = sigma_mag_func.vector().get_local()
    mu_array = mu.vector().get_local()
    lmbda_array = lmbda.vector().get_local()
    rho_gfap_arr = rho_gfap.vector().get_local()
    rho_ecm_arr = rho_ecm.vector().get_local()
    J_g_array = J_g_func.vector().get_local()
    Ic_e_array = Ic_e_func.vector().get_local()
    theta_array = theta_func.vector().get_local()
    
    # Mechanobiological Stiffening (ODE integration)
    Phi_arr = np.maximum(0.0, (sigma_mag_array - sigma_homeo) / sigma_homeo)
    rho_gfap_arr += dt_growth * (k_gfap * Phi_arr - d_gfap * (rho_gfap_arr - 1.0))
    rho_ecm_arr += dt_growth * (k_ecm * Phi_arr - d_ecm * (rho_ecm_arr - 1.0))
    
    # Rule of Mixtures: Intrinsic baseline stiffness is augmented by scar tissue synthesis
    mu_array = mu_val + c1 * (rho_gfap_arr - 1.0) + c2 * (rho_ecm_arr - 1.0)
    lmbda_array = lmbda_val * (mu_array / mu_val) # Scale bulk modulus proportionally
    
    p_cap_MPa = 25.0 * 0.0001333
    Pt_array = p_cap_MPa + p_tissue_array # Add because p_tissue_array (min principal stress) is negative in compression
    
    P_crit = 14 * 0.0001333  # Increased capillary fragility threshold
    k = 2.0 / 0.0001333
    
    # Calculate local blood flow reduction
    Q_local_array = 0.5 + 0.5 * (1.0 / (1.0 + np.exp(-k * (Pt_array - P_crit))))
    
    # Compute local atrophy trigger (only shrink if flow < 0.74)
    atrophy_array = np.where(Q_local_array < 0.74, 0.74 - Q_local_array, 0.0)
    
    # --- 3. Update FEniCS Fields ---
    mu.vector().set_local(mu_array)
    mu.vector().apply("insert")
    lmbda.vector().set_local(lmbda_array)
    lmbda.vector().apply("insert")
    rho_gfap.vector().set_local(rho_gfap_arr)
    rho_gfap.vector().apply("insert")
    rho_ecm.vector().set_local(rho_ecm_arr)
    rho_ecm.vector().apply("insert")

    atrophy_field.vector().set_local(atrophy_array)
    atrophy_field.vector().apply("insert")
    
    # --- Atrophy/Shrinkage increment ---
    theta_array -= growth_rate * atrophy_array * dt_growth * theta_array
    theta_array = np.clip(theta_array, 0.1, 1.0)
    theta_func.vector().set_local(theta_array)
    theta_func.vector().apply("insert")

    # --- Solve ---
    # The solver 're-reads' the updated theta values inside the symbolic F_res_growth form
    solve(F_res_growth == 0, u, bcs=[bc_sym_axis, bc_bottom_roller, bc_outer_roller], J=J_res_growth, solver_parameters=solver_params)

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
    xdmf_growth.write(p_func, float(gstep))
    xdmf_growth.write(hoop_func, float(gstep))
    xdmf_growth.write(sig1_func, float(gstep))
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
Fe_un = F_current_un * inv(Fg_ufl_un)
C_un = Fe_un.T * Fe_un
Ic_un = tr(C_un)
J_un = det(Fe_un)
psi_un = (mu / 2) * (Ic_un - 3) - mu * ln(J_un) + (lmbda / 2) * (ln(J_un))**2
Pi_total_unload = psi_un * r_fac * dx_c + load_factor * iop_target_MPa * dot(n, u) * r_fac * ds_measure(pressure_id)
F_res_unload = derivative(Pi_total_unload, u, v)
J_res_unload = derivative(F_res_unload, u, du)

# UFL for exact Cauchy stress during unloading
F_var_U = variable(get_F3D(u))
Fg_ufl_U = sqrt(theta_func) * I3 + (1.0 - sqrt(theta_func)) * outer(n0, n0)
Fe_var_U = F_var_U * inv(Fg_ufl_U)
C_e_U = Fe_var_U.T * Fe_var_U
Ic_e_U = tr(C_e_U)
J_e_U = det(Fe_var_U)
J_tot_U = det(F_var_U)
Psi_U = (mu / 2) * (Ic_e_U - 3) - mu * ln(J_e_U) + (lmbda / 2) * (ln(J_e_U))**2
P_exact_U = diff(Psi_U, F_var_U)
sigma_exact_U = (1.0 / J_tot_U) * P_exact_U * F_var_U.T
sigma_dev_U = sigma_exact_U - (1.0 / 3.0) * tr(sigma_exact_U) * I3
vm_expr_U = sqrt((3.0 / 2.0) * inner(sigma_dev_U, sigma_dev_U))

# Calculate Minimum Principal Stress (Max Compressive)
sig_center_U = (sigma_exact_U[0, 0] + sigma_exact_U[2, 2]) / 2.0
sig_radius_U = sqrt(((sigma_exact_U[0, 0] - sigma_exact_U[2, 2]) / 2.0)**2 + sigma_exact_U[0, 2]**2)
sig_3_U = conditional(lt(sig_center_U - sig_radius_U, sigma_exact_U[1, 1]), sig_center_U - sig_radius_U, sigma_exact_U[1, 1])
p_expr_U = sig_3_U
hoop_expr_U = sigma_exact_U[1, 1]
sig_1_U = conditional(gt(sig_center_U + sig_radius_U, sigma_exact_U[1, 1]), sig_center_U + sig_radius_U, sigma_exact_U[1, 1])

print("=== Unloading Phase ===")
for step in range(1, unload_steps + 1):
    # Gradually reduce load and scleral expansion
    factor = 1.0 * (1 - step / unload_steps)
    load_factor.assign(factor) # Scale from 1.0 down to 0.0
    u_outer_rad.assign(0.05 * x_max * factor)

    # --- Solve equilibrium ---
    solve(F_res_unload == 0, u, bcs=[bc_sym_axis, bc_bottom_roller, bc_outer_roller], J=J_res_unload, solver_parameters=solver_params)
    
    vm_func.assign(project(vm_expr_U, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    p_func.assign(project(p_expr_U, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    hoop_func.assign(project(hoop_expr_U, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    sig1_func.assign(project(sig_1_U, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    
    xdmf_unload.write(u, float(step))
    xdmf_unload.write(vm_func, float(step))
    xdmf_unload.write(p_func, float(step))
    xdmf_unload.write(hoop_func, float(step))
    xdmf_unload.write(sig1_func, float(step))

xdmf_unload.close()
print("✅ Unloading phase finished and saved.")
