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

V = VectorFunctionSpace(mesh, "CG", 2)  # 2 means quadratic

# ======================================================
# 3. Boundary conditions and material parameters
# ======================================================
iop_target_MPa = 40 * 0.0001333  
csfp_target_MPa = 0.00 * 0.0001333 

# Lamina cribrosa compliance derivation (Winkler elastic foundation)
E_LC_MPa = 0.3  
t_LC_mm = 0.3   
k_LC = Constant(E_LC_MPa / t_LC_mm) 

zero = Constant((0.0, 0.0))

axis_id = ID_SYMMETRY_AXIS
bc_sym_axis = DirichletBC(V.sub(0), Constant(0.0), boundaries, axis_id)

bottom_robin_id = ID_LC_FLOOR

u_outer_rad = Constant(0.0) 

bc_sclera_wall_x = DirichletBC(V.sub(0), u_outer_rad, boundaries, ID_SCLERA_WALL)
bc_sclera_wall_y = DirichletBC(V.sub(1), Constant(0.0), boundaries, ID_SCLERA_WALL)
bc_retina_roller_y = DirichletBC(V.sub(1), Constant(0.0), boundaries, ID_RETINA_BOUNDARY)

all_bcs = [bc_sym_axis, bc_sclera_wall_x, bc_sclera_wall_y, bc_retina_roller_y]

nu = 0.499                                                                          
mu_val = 56.0 / 1000                                                                
E_val = 2 * mu_val * (1 + nu)                                                       
lmbda_val = (E_val * nu) / ((1 + nu) * (1 - 2 * nu)) 

V_DG0_scalar = FunctionSpace(mesh, "DG", 0)
V_DG0_tensor = TensorFunctionSpace(mesh, "DG", 0, shape=(3, 3))

# --- PURE CMT CONSTITUENT INTRINSIC PROPERTIES ---
mu_axon_pure  = Constant(mu_val * 0.2)   # Soft cellular neural channels
mu_astro_pure = Constant(mu_val * 0.5)   # Glial structural baseline processes
mu_ecm_pure   = Constant(mu_val * 3.0)   # Stiff structural matrix scaffolding

# --- UPGRADE: ADD ANISOTROPIC FIBER STIFFNESS PARAMETERS ---
mu_axon_fiber  = Constant(mu_val * 0.4)   # Directional resistance of neural structural channels
mu_astro_fiber = Constant(mu_val * 0.8)   # Directional resistance of wrapping glial tubes

u = Function(V)
u.rename("Displacement", "Displacement")
du = TrialFunction(V)
v = TestFunction(V)

I3 = Identity(3)

# Analytical smooth maximum to replace sharp conditional jumps and ensure Newton convergence
def smooth_max(a, b, eps=1e-6):
    return 0.5 * (a + b + sqrt((a - b)**2 + eps))

def get_F3D(u):
    r = SpatialCoordinate(mesh)[0]
    r_safe = r + 1e-14 
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

dx_c = Measure("dx", domain=mesh, metadata={"quadrature_degree": 4})
r_coord = SpatialCoordinate(mesh)[0]
r_fac = 2 * np.pi * r_coord

# ======================================================
# 4. Loading Phase (step load)
# ======================================================
load_steps = 100
pressure_id = ID_CUP_SURFACE_IOP 
load_factor = Constant(0.0)
ds_measure = Measure("ds", domain=mesh, subdomain_data=boundaries)
n = FacetNormal(mesh)

F = get_F3D(u)
C = F.T * F
Ic = tr(C)
J = det(F)

psi = (Constant(mu_val) / 2) * (Ic - 3) - Constant(mu_val) * ln(J) + (Constant(lmbda_val) / 2) * (ln(J))**2   

Pi_total = (psi * r_fac * dx_c 
            + load_factor * iop_target_MPa * dot(n, u) * r_fac * ds_measure(pressure_id)
            + load_factor * csfp_target_MPa * dot(n, u) * r_fac * ds_measure(bottom_robin_id)
            + (k_LC / 2.0) * (u[1]**2) * r_fac * ds_measure(bottom_robin_id))
F_res = derivative(Pi_total, u, v)
J_res = derivative(F_res, u, du)

solver_params = {"newton_solver": {"relative_tolerance": 1e-8, 
                                   "absolute_tolerance": 1e-10, 
                                   "maximum_iterations": 25, 
                                   "linear_solver": "mumps"}}

xdmf_load = XDMFFile(mesh.mpi_comm(), "u_loading_phase.xdmf")
xdmf_load.parameters["flush_output"] = True
xdmf_load.parameters["functions_share_mesh"] = True

F_var_L = variable(get_F3D(u))
C_L = F_var_L.T * F_var_L
Ic_L = tr(C_L)
J_L = det(F_var_L)
Psi_L = (Constant(mu_val) / 2) * (Ic_L - 3) - Constant(mu_val) * ln(J_L) + (Constant(lmbda_val) / 2) * (ln(J_L))**2
P_L = diff(Psi_L, F_var_L)
sigma_L = (1.0 / J_L) * P_L * F_var_L.T
sigma_dev_L = sigma_L - (1.0 / 3.0) * tr(sigma_L) * I3
vm_expr_L = sqrt((3.0 / 2.0) * inner(sigma_dev_L, sigma_dev_L))

sig_center_L = (sigma_L[0, 0] + sigma_L[2, 2]) / 2.0
sig_radius_L = sqrt(((sigma_L[0, 0] - sigma_L[2, 2]) / 2.0)**2 + sigma_L[0, 2]**2)
sig_3_L = conditional(lt(sig_center_L - sig_radius_L, sigma_L[1, 1]), sig_center_L - sig_radius_L, sigma_L[1, 1])
p_expr_L = sig_3_L
hoop_expr_L = sigma_L[1, 1]

E_L = 0.5 * (C_L - I3)
E_center_L = (E_L[0, 0] + E_L[2, 2]) / 2.0
E_radius_L = sqrt(((E_L[0, 0] - E_L[2, 2]) / 2.0)**2 + E_L[0, 2]**2)
E_1_L = conditional(gt(E_center_L + E_radius_L, E_L[1, 1]), E_center_L + E_radius_L, E_L[1, 1])
E_3_L = conditional(lt(E_center_L - E_radius_L, E_L[1, 1]), E_center_L - E_radius_L, E_L[1, 1])

R_eye = 12.0 
t_sclera = 1.0 
E_sclera = 3.0 

sigma_0_max = (iop_target_MPa * R_eye) / (2.0 * t_sclera)
sigma_theta_max = 2.0 * sigma_0_max
epsilon_theta_max = (1.0 / E_sclera) * sigma_theta_max
canal_radius = 0.8 
Delta_r_max = canal_radius * epsilon_theta_max

print("=== Loading Phase ===")
for step in range(1, load_steps + 1):
    factor = step / load_steps
    load_factor.assign(factor) 
    
    u_outer_rad.assign(Delta_r_max * factor)
    solve(F_res == 0, u, bcs=all_bcs, J=J_res, solver_parameters=solver_params)
    
    vm_func.assign(project(vm_expr_L, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    E3_func.assign(project(E_3_L, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    hoop_func.assign(project(hoop_expr_L, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    doc = project(E_1_L, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2})
    E1_func.assign(doc)
    
    xdmf_load.write(u, float(step))
    xdmf_load.write(vm_func, float(step))
    xdmf_load.write(E3_func, float(step))
    xdmf_load.write(hoop_func, float(step))
    xdmf_load.write(E1_func, float(step))
xdmf_load.close()
print("✅ Loading phase finished and saved.")

F_deposit = project(get_F3D(u), V_DG0_tensor, form_compiler_parameters={"quadrature_degree": 2})

# ======================================================
# 5. Remodeling Phase (maintain load)
# ======================================================
growth_steps = 500
dt_growth = 0.1
growth_rate = 0.05

# Initialize separate physical volume fractions fields
theta_axon  = Function(V_DG0_scalar)
theta_axon.vector()[:] = 0.45

theta_astro = Function(V_DG0_scalar)
theta_astro.vector()[:] = 0.30

theta_ecm   = Function(V_DG0_scalar)
theta_ecm.vector()[:] = 0.25

theta_total = Function(V_DG0_scalar)
theta_total.rename("Theta_Total", "Current Volume Fraction")
theta_total.vector()[:] = 1.0

# Fields for tracking multi-configuration snapshot adjustments
theta_ecm_old = Function(V_DG0_scalar)
theta_ecm_old.vector()[:] = 0.25
delta_theta_ecm = Function(V_DG0_scalar)
delta_theta_ecm.vector()[:] = 0.0

# --- Solve Laplace PDE for Axon Orientation ---
V_phi = FunctionSpace(mesh, "CG", 2)
phi = TrialFunction(V_phi)
v_phi = TestFunction(V_phi)

a_phi = inner(grad(phi), grad(v_phi)) * r_fac * dx_c
L_phi = Constant(0.0) * v_phi * r_fac * dx_c

bc_source = DirichletBC(V_phi, Constant(1.0), boundaries, ID_RETINA_BOUNDARY)
bc_sink = DirichletBC(V_phi, Constant(0.0), boundaries, ID_LC_FLOOR)

phi_sol = Function(V_phi)
solve(a_phi == L_phi, phi_sol, [bc_source, bc_sink])

grad_phi = grad(phi_sol)
grad_norm = sqrt(inner(grad_phi, grad_phi)) + 1e-14 
n0_2d = -grad_phi / grad_norm
n0 = as_vector((n0_2d[0], 0.0, n0_2d[1]))

V_DG0_vec = VectorFunctionSpace(mesh, "DG", 0)
n0_func = project(n0_2d, V_DG0_vec, form_compiler_parameters={"quadrature_degree": 2})
n0_func.rename("Axon_Orientation", "Reference Axon Direction")

with XDMFFile(mesh.mpi_comm(), "axon_orientation.xdmf") as xdmf_n0:
    xdmf_n0.write(n0_func)

load_factor.assign(1.0) 
u_outer_rad.assign(Delta_r_max)

xdmf_growth = XDMFFile(mesh.mpi_comm(), "u_growth_phase.xdmf")
xdmf_growth.parameters["flush_output"] = True
xdmf_growth.parameters["functions_share_mesh"] = True

xdmf_theta = XDMFFile(mesh.mpi_comm(), "theta_growth_phase.xdmf")
xdmf_theta.parameters["flush_output"] = True
xdmf_theta.parameters["functions_share_mesh"] = True

# ======================================================
# UPGRADED: ANISOTROPIC CMT VARIATIONAL ASSEMBLY
# ======================================================
F_current_sym = get_F3D(u)
C_sym = F_current_sym.T * F_current_sym
Ic_sym = tr(C_sym)
J_sym = det(F_current_sym)

Fe_scar_sym = F_current_sym * inv(F_deposit)
C_scar_sym = Fe_scar_sym.T * Fe_scar_sym
Ic_scar_sym = tr(C_scar_sym)
J_scar_sym = det(Fe_scar_sym)

# --- Add Anisotropic Structural Invariants to the Variational Equations ---
I4_sym = dot(n0, C_sym * n0)
I4_sym_gated = smooth_max(I4_sym, 1.0, eps=1e-6)  # Tension-only structural safety gate

# Mechanical blocks receive the + (mu_fiber/2)*(I4 - 1)**2 directional resistance additions
psi_shear_axon  = theta_axon  * ((mu_axon_pure / 2)  * (Ic_sym - 3 - 2 * ln(J_sym)) + (mu_axon_fiber / 2) * (I4_sym_gated - 1)**2)
psi_shear_astro = theta_astro * ((mu_astro_pure / 2) * (Ic_sym - 3 - 2 * ln(J_sym)) + (mu_astro_fiber / 2) * (I4_sym_gated - 1)**2)
psi_shear_ecm   = Constant(0.25) * (mu_ecm_pure / 2)   * (Ic_sym - 3 - 2 * ln(J_sym))
psi_shear_scar  = (theta_ecm - Constant(0.25)) * (mu_ecm_pure / 2) * (Ic_scar_sym - 3 - 2 * ln(J_scar_sym))
psi_volumetric  = (Constant(lmbda_val) / 2) * (ln(J_sym / theta_total))**2

psi_growth = psi_shear_axon + psi_shear_astro + psi_shear_ecm + psi_shear_scar + psi_volumetric

Pi_total_growth = (psi_growth * r_fac * dx_c 
                    + load_factor * iop_target_MPa * dot(n, u) * r_fac * ds_measure(pressure_id)
                    + load_factor * csfp_target_MPa * dot(n, u) * r_fac * ds_measure(bottom_robin_id)
                    + (k_LC / 2.0) * (u[1]**2) * r_fac * ds_measure(bottom_robin_id))
F_res_growth = derivative(Pi_total_growth, u, v)
J_res_growth = derivative(F_res_growth, u, du)

print("=== Remodeling Phase (Atrophy/Shrinkage + Stiffening) ===")

log_file = open("remodeling_log.txt", "w")
log_file.write("Step,Atrophy Rate,Avg Apparent Modulus,Max Apparent Modulus\n")

# Biological ODE Parameters
strain_homeo = 1.8e-2  
k_astro = 0.0010      
d_astro = 0.20        
k_ecm_dep = 0.0010    
d_ecm = 0.20         

theta_min = 0.3      
gamma_atrophy = 2.0  

atrophy_field = Function(V_DG0_scalar)
atrophy_field.rename("Atrophy_Signal", "Atrophy_Signal")
shrinkage_field = Function(V_DG0_scalar)
shrinkage_field.rename("Shrinkage", "Shrinkage")

# ======================================================
# UPGRADED: ANISOTROPIC CAUCHY STRESS TENSOR EXTRACTION
# ======================================================
F_var = variable(get_F3D(u))
C_var = F_var.T * F_var
Ic_var = tr(C_var)
J_var = det(F_var)

Fe_scar_var = F_var * inv(F_deposit)
C_scar_var = Fe_scar_var.T * Fe_scar_var
Ic_scar_var = tr(C_scar_var)
J_scar_var = det(Fe_scar_var)

# --- Add Anisotropic Structural Invariants to the Post-Processing Stress Matrix ---
I4_var = dot(n0, C_var * n0)
I4_var_gated = smooth_max(I4_var, 1.0, eps=1e-6)

psi_axon_v  = theta_axon  * ((mu_axon_pure / 2)  * (Ic_var - 3 - 2 * ln(J_var)) + (mu_axon_fiber / 2) * (I4_var_gated - 1)**2)
psi_astro_v = theta_astro * ((mu_astro_pure / 2) * (Ic_var - 3 - 2 * ln(J_var)) + (mu_astro_fiber / 2) * (I4_var_gated - 1)**2)
psi_ecm_v   = Constant(0.25) * (mu_ecm_pure / 2)   * (Ic_var - 3 - 2 * ln(J_var))
psi_scar_v  = (theta_ecm - Constant(0.25)) * (mu_ecm_pure / 2) * (Ic_scar_var - 3 - 2 * ln(J_scar_var))
psi_vol_v   = (Constant(lmbda_val) / 2) * (ln(J_var / theta_total))**2

Psi_tot_var = psi_axon_v + psi_astro_v + psi_ecm_v + psi_scar_v + psi_vol_v
P_exact = diff(Psi_tot_var, F_var)
J_tot = det(F_var)
sigma_exact = (1.0 / J_tot) * P_exact * F_var.T

# Extract Directional Strains via Local Basis Transformations
E_neural_var = 0.5 * (C_var - I3)

m0 = as_vector((-n0[2], 0.0, n0[0]))          
e_theta = as_vector((0.0, 1.0, 0.0))          

E_parallel = dot(n0, E_neural_var * n0)        
E_m        = dot(m0, E_neural_var * m0)        
E_theta    = E_neural_var[1, 1]                

E_perp_min = conditional(lt(E_m, E_theta), E_m, E_theta)  

sigma_dev_exact = sigma_exact - (1.0 / 3.0) * tr(sigma_exact) * I3
vm_expr_G = sqrt((3.0 / 2.0) * inner(sigma_dev_exact, sigma_dev_exact))

J_g_expr = Constant(1.0)
strain_mag_expr = E_parallel  

for gstep in range(1, growth_steps + 1):
    # --- Projections ---
    E3_func.assign(project(E_perp_min, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    vm_func.assign(project(vm_expr_G, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    hoop_func.assign(project(sigma_exact[1, 1], V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    E1_func.assign(project(E_parallel, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    
    strain_mag_func = project(strain_mag_expr, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2})
    J_g_func = project(J_g_expr, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2})
    Ic_e_func = project(Ic_var, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2})
    
    theta_ecm_old.assign(theta_ecm)
    
    # --- 2. Fast NumPy Biological Math ---
    E3_array = E3_func.vector().get_local()
    strain_mag_array = strain_mag_func.vector().get_local()
    Ic_e_array = Ic_e_func.vector().get_local()
    
    theta_axon_arr  = theta_axon.vector().get_local()
    theta_astro_arr = theta_astro.vector().get_local()
    theta_ecm_arr   = theta_ecm.vector().get_local()

    Phi_arr = np.maximum(0.0, (strain_mag_array - strain_homeo) / strain_homeo)
    
    theta_astro_arr += dt_growth * (k_astro * Phi_arr - d_astro * (theta_astro_arr - 0.30))
    theta_ecm_arr   += dt_growth * (k_ecm_dep * Phi_arr - d_ecm * (theta_ecm_arr - 0.25))
    
    delta_theta_ecm_arr = theta_ecm_arr - theta_ecm_old.vector().get_local()
    
    lambda3_sq = np.clip(2.0 * E3_array + 1.0, 0.0, None)
    Q_local_array = lambda3_sq**2
    atrophy_array = np.where(Q_local_array < 0.876096, 0.876096 - Q_local_array, 0.0)
    
    axon_normalized = theta_axon_arr / 0.45
    k_theta = growth_rate * ((axon_normalized - theta_min) / (1.0 - theta_min))**gamma_atrophy
    k_theta = np.maximum(k_theta, 0.0)
    
    theta_axon_arr -= k_theta * atrophy_array * dt_growth * theta_axon_arr
    theta_axon_arr = np.clip(theta_axon_arr, theta_min * 0.45, 0.45)
    
    # --- 3. Update FEniCS Fields ---
    theta_axon.vector().set_local(theta_axon_arr)
    theta_axon.vector().apply("insert")
    theta_astro.vector().set_local(theta_astro_arr)
    theta_astro.vector().apply("insert")
    theta_ecm.vector().set_local(theta_ecm_arr)
    theta_ecm.vector().apply("insert")
    delta_theta_ecm.vector().set_local(delta_theta_ecm_arr)
    delta_theta_ecm.vector().apply("insert")
    
    atrophy_field.vector().set_local(atrophy_array)
    atrophy_field.vector().apply("insert")

    theta_total_arr = theta_axon_arr + theta_astro_arr + theta_ecm_arr
    theta_total.vector().set_local(theta_total_arr)
    theta_total.vector().apply("insert")

    shrinkage_array = 1.0 - theta_total_arr
    shrinkage_field.vector().set_local(shrinkage_array)
    shrinkage_field.vector().apply("insert")

    theta_ecm_new_mass = theta_ecm - Constant(0.25)
    theta_ecm_new_mass_safe = conditional(lt(theta_ecm_new_mass, 1e-12), 1e-12, theta_ecm_new_mass)
    
    F_deposit_expr = conditional(gt(delta_theta_ecm, 0.0),
                                 ((theta_ecm_old - Constant(0.25)) * F_deposit + delta_theta_ecm * get_F3D(u)) / theta_ecm_new_mass_safe,
                                 F_deposit)
    F_deposit.assign(project(F_deposit_expr, V_DG0_tensor, form_compiler_parameters={"quadrature_degree": 2}))

    solve(F_res_growth == 0, u, bcs=all_bcs, J=J_res_growth, solver_parameters=solver_params)

    if gstep % 10 == 0:
        mu_apparent = (theta_axon_arr * mu_val * 0.2 + theta_astro_arr * mu_val * 0.5 + 0.25 * mu_val * 3.0 + (theta_ecm_arr - 0.25) * mu_val * 3.0)
        avg_atrophy = np.mean(atrophy_array)
        avg_mu = np.mean(mu_apparent)
        max_mu = np.max(mu_apparent)
        print(f"[Remodeling Step {gstep}] Atrophy Rate = {avg_atrophy:.4e}")
        print(f"                     Avg Apparent mu = {avg_mu:.4f} (Max: {max_mu:.4f})")
        print(f"                     Max I_c (Distortion) = {np.max(Ic_e_array):.4f}")
        
        log_line = f"{gstep},{avg_atrophy:.4e},{avg_mu:.4f},{max_mu:.4f}\n"
        log_file.write(log_line)

    xdmf_growth.write(u, float(gstep))
    xdmf_growth.write(vm_func, float(gstep))
    xdmf_growth.write(E3_func, float(gstep))
    xdmf_growth.write(hoop_func, float(gstep))
    xdmf_growth.write(E1_func, float(gstep))
    xdmf_growth.write(atrophy_field, float(gstep))
    xdmf_growth.write(shrinkage_field, float(gstep))
    if gstep % 20 == 0:
        xdmf_theta.write(theta_total, float(gstep))

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

F_current_un = get_F3D(u)
C_un = F_current_un.T * F_current_un
Ic_un = tr(C_un)
J_un = det(F_current_un)

Fe_scar_un = F_current_un * inv(F_deposit)
C_scar_un = Fe_scar_un.T * Fe_scar_un
Ic_scar_un = tr(C_scar_un)
J_scar_un = det(Fe_scar_un)

# --- Add Anisotropic Structural Invariants to Unloading Variational Assembly ---
I4_un = dot(n0, C_un * n0)
I4_un_gated = smooth_max(I4_un, 1.0, eps=1e-6)

psi_shear_axon_un  = theta_axon  * ((mu_axon_pure / 2)  * (Ic_un - 3 - 2 * ln(J_un)) + (mu_axon_fiber / 2) * (I4_un_gated - 1)**2)
psi_shear_astro_un = theta_astro * ((mu_astro_pure / 2) * (Ic_un - 3 - 2 * ln(J_un)) + (mu_astro_fiber / 2) * (I4_un_gated - 1)**2)
psi_shear_ecm_un   = Constant(0.25) * (mu_ecm_pure / 2)   * (Ic_un - 3 - 2 * ln(J_un))
psi_shear_scar_un  = (theta_ecm - Constant(0.25)) * (mu_ecm_pure / 2) * (Ic_scar_un - 3 - 2 * ln(J_scar_un))
psi_volumetric_un  = (Constant(lmbda_val) / 2) * (ln(J_un / theta_total))**2

psi_un = psi_shear_axon_un + psi_shear_astro_un + psi_shear_ecm_un + psi_shear_scar_un + psi_volumetric_un

Pi_total_unload = (psi_un * r_fac * dx_c 
                    + load_factor * iop_target_MPa * dot(n, u) * r_fac * ds_measure(pressure_id)
                    + load_factor * csfp_target_MPa * dot(n, u) * r_fac * ds_measure(bottom_robin_id)
                    + (k_LC / 2.0) * (u[1]**2) * r_fac * ds_measure(bottom_robin_id))
F_res_unload = derivative(Pi_total_unload, u, v)
J_res_unload = derivative(F_res_unload, u, du)

F_var_U = variable(get_F3D(u))
C_U = F_var_U.T * F_var_U
Ic_U = tr(C_U)
J_U = det(F_var_U)

Fe_scar_U = F_var_U * inv(F_deposit)
C_scar_U = Fe_scar_U.T * Fe_scar_U
Ic_scar_U = tr(C_scar_U)
J_scar_U = det(Fe_scar_U)

# --- Add Anisotropic Structural Invariants to Unloading Post-Processing Stress Matrix ---
I4_U = dot(n0, C_U * n0)
I4_U_gated = smooth_max(I4_U, 1.0, eps=1e-6)

psi_axon_U  = theta_axon  * ((mu_axon_pure / 2)  * (Ic_U - 3 - 2 * ln(J_U)) + (mu_axon_fiber / 2) * (I4_U_gated - 1)**2)
psi_astro_U = theta_astro * ((mu_astro_pure / 2) * (Ic_U - 3 - 2 * ln(J_U)) + (mu_astro_fiber / 2) * (I4_U_gated - 1)**2)
psi_ecm_U   = Constant(0.25) * (mu_ecm_pure / 2)   * (Ic_var - 3 - 2 * ln(J_var)) # Anchored safely to static load invariants
psi_scar_U  = (theta_ecm - Constant(0.25)) * (mu_ecm_pure / 2) * (Ic_scar_U - 3 - 2 * ln(J_scar_U))
psi_vol_U   = (Constant(lmbda_val) / 2) * (ln(J_U / theta_total))**2

Psi_U = psi_axon_U + psi_astro_U + psi_ecm_U + psi_scar_U + psi_vol_U
P_exact_U = diff(Psi_U, F_var_U)
J_tot_U = det(F_var_U)
sigma_exact_U = (1.0 / J_tot_U) * P_exact_U * F_var_U.T

sigma_dev_U = sigma_exact_U - (1.0 / 3.0) * tr(sigma_exact_U) * I3
vm_expr_U = sqrt((3.0 / 2.0) * inner(sigma_dev_U, sigma_dev_U))
hoop_expr_U = sigma_exact_U[1, 1]

E_neural_U = 0.5 * (C_U - I3)
E_parallel_U = dot(n0, E_neural_U * n0)
E_m_U        = dot(m0, E_neural_U * m0)
E_theta_U    = E_neural_U[1, 1]
E_perp_min_U = conditional(lt(E_m_U, E_theta_U), E_m_U, E_theta_U)

print("=== Unloading Phase ===")
for step in range(1, unload_steps + 1):
    factor = 1.0 * (1 - step / unload_steps)
    load_factor.assign(factor) 
    u_outer_rad.assign(Delta_r_max * factor)

    solve(F_res_unload == 0, u, bcs=all_bcs, J=J_res_unload, solver_parameters=solver_params)
    
    vm_func.assign(project(vm_expr_U, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    E3_func.assign(project(E_perp_min_U, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    hoop_func.assign(project(hoop_expr_U, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    E1_func.assign(project(E_parallel_U, V_DG0_scalar, form_compiler_parameters={"quadrature_degree": 2}))
    
    xdmf_unload.write(u, float(step))
    xdmf_unload.write(vm_func, float(step))
    xdmf_unload.write(E3_func, float(step))
    xdmf_unload.write(hoop_func, float(step))
    xdmf_unload.write(E1_func, float(step))

xdmf_unload.close()
print("✅ Unloading phase finished and saved.")