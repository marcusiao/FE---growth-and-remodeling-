import os
from dolfin import *
import numpy as np

# ==========================================================
# 1. Setup Path and Load 3D Mesh, Subdomains & Boundaries
# ==========================================================
os.chdir(os.path.dirname(os.path.realpath(__file__)))

file_name = "3D_ONH_corase"

print(f"Loading {file_name} 3D mesh, subdomains, and boundary facets...")

mesh = Mesh()
with XDMFFile(f"{file_name}_domain.xdmf") as infile:
    infile.read(mesh)

# Subdomains (1: Choroid, 2: LC, 3: PT, 4: Sclera)
mvc_domain = MeshValueCollection("size_t", mesh, mesh.topology().dim())
with XDMFFile(f"{file_name}_domain.xdmf") as infile:
    infile.read(mvc_domain, "name_to_read")
subdomains = MeshFunction("size_t", mesh, mvc_domain)

# Boundary Facets (1: SYMMETRY, 2: IOP, 3: FIXED, 4: OUT)
mvc_facets = MeshValueCollection("size_t", mesh, mesh.topology().dim() - 1)
with XDMFFile(f"{file_name}_facet_region.xdmf") as infile:
    infile.read(mvc_facets, "name_to_read")
facet_markers = MeshFunction("size_t", mesh, mvc_facets)

ds_measure = Measure("ds", domain=mesh, subdomain_data=facet_markers)

# ==========================================================
# 2. Function Spaces & Growth Field Initialization
# ==========================================================
V = VectorFunctionSpace(mesh, "CG", 1)
V_DG0 = FunctionSpace(mesh, "DG", 0)
V_tensor = TensorFunctionSpace(mesh, "DG", 0, shape=(3, 3))

u = Function(V)
u.rename("Displacement", "u")
du = TrialFunction(V)
v = TestFunction(V)

# Growth Scalar Field (theta = 1.0 is healthy reference state)
theta_func = Function(V_DG0)
theta_func.rename("Theta_Growth", "theta")
theta_func.vector()[:] = 1.0

# ==========================================================
# 3. Material Properties Mapping (DG0 Fields)
# ==========================================================
material_props = {
    1: {"E": 0.2, "nu": 0.49},  # Choroid
    2: {"E": 0.3, "nu": 0.49},  # Lamina Cribrosa (LC)
    3: {"E": 0.1, "nu": 0.49},  # Prelaminar Tissue (PT)
    4: {"E": 3.0, "nu": 0.49},  # Sclera
}

mu_field = Function(V_DG0)
lmbda_field = Function(V_DG0)

sub_arr = subdomains.array()
mu_arr = np.zeros(mesh.num_cells())
lmbda_arr = np.zeros(mesh.num_cells())

for tag, props in material_props.items():
    E_val = props["E"]
    nu_val = props["nu"]
    mu_val = E_val / (2.0 * (1.0 + nu_val))
    lmbda_val = (E_val * nu_val) / ((1.0 + nu_val) * (1.0 - 2.0 * nu_val))

    mask = sub_arr == tag
    mu_arr[mask] = mu_val
    lmbda_arr[mask] = lmbda_val

dof_map = V_DG0.dofmap()
for cell_idx in range(mesh.num_cells()):
    dof = dof_map.cell_dofs(cell_idx)[0]
    mu_field.vector()[dof] = mu_arr[cell_idx]
    lmbda_field.vector()[dof] = lmbda_arr[cell_idx]

# ==========================================================
# 4. Kinematics (F = Fe * Fg) & Hyperelastic Potential
# ==========================================================
I3 = Identity(3)
F = I3 + grad(u)
Fe = F / theta_func
Ce = Fe.T * Fe
Ic_e = tr(Ce)
Je = det(Fe)

psi = (mu_field / 2.0) * (Ic_e - 3.0) - mu_field * ln(Je) + (lmbda_field / 2.0) * (ln(Je)) ** 2

# ==========================================================
# 5. Realistic Boundary Conditions & Variational Formulation
# ==========================================================
# 1. FIXED (Tag 3) -> Constrain Y-displacement only (roller on X-Z plane)
bc_fixed_y = DirichletBC(V.sub(1), Constant(0.0), facet_markers, 3)

# 2. SYMMETRY (Tag 1) -> Fix Z displacement only (allows free in-plane X-Y deformation)
bc_sym_z = DirichletBC(V.sub(2), Constant(0.0), facet_markers, 1)

# 3. OUT Boundary (Tag 4) -> Dynamic radial expansion in X-Z plane
# Link scleral expansion to IOP via analytical model (adapted from 2D)
R_eye = 12.0      # mm, radius of the eyeball
t_sclera = 1.0    # mm, thickness of the sclera
E_sclera = 3.0    # MPa, Young's modulus of sclera
canal_radius = 0.8 # mm, inner radius of the scleral canal

iop_target_MPa = 0.005  # Target IOP (~37.5 mmHg)

# Calculate max expansion based on target IOP
sigma_0_max = (iop_target_MPa * R_eye) / (2.0 * t_sclera) # Laplace's Law
sigma_theta_max = 2.0 * sigma_0_max # Stress concentration
epsilon_theta_max = (1.0 / E_sclera) * sigma_theta_max # Hooke's Law
Delta_r_max = canal_radius * epsilon_theta_max

disp_factor = Constant(0.0) # Use a Constant for robust updates

rad_disp_expr = Expression(
    ("factor * Delta_r * x[0] / (sqrt(x[0]*x[0] + x[2]*x[2]) + 1e-12)",
     "0.0",
     "factor * Delta_r * x[2] / (sqrt(x[0]*x[0] + x[2]*x[2]) + 1e-12)"),
    factor=disp_factor,
    Delta_r=Delta_r_max,
    degree=1)

bc_out_radial = DirichletBC(V, rad_disp_expr, facet_markers, 4)

bcs = [bc_fixed_y, bc_sym_z, bc_out_radial]

# 4. IOP Normal Pressure (Tag 2)
load_factor = Constant(0.0)
n_facet = FacetNormal(mesh)

# Potential Energy: Strain Energy + Normal Pressure Work
Pi = psi * dx + load_factor * iop_target_MPa * dot(n_facet, u) * ds_measure(2)
F_res = derivative(Pi, u, v)
J_res = derivative(F_res, u, du)

# Solver Settings
solver_params = {
    "newton_solver": {
        "linear_solver": "mumps",
        "relative_tolerance": 1e-6,
        "absolute_tolerance": 1e-8,
        "maximum_iterations": 25,
    }
}

# ==========================================================
# 6. Post-Processing Quantities
# ==========================================================
E_gl_expr = 0.5 * (Ce - I3)
E1_func = Function(V_DG0)
E1_func.rename("Max_Principal_Strain", "E1")

vm_func = Function(V_DG0)
vm_func.rename("von_Mises_Stress", "vm")

P_expr = diff(psi, variable(F))
sigma_expr = (1.0 / det(F)) * P_expr * F.T
sigma_dev = sigma_expr - (1.0 / 3.0) * tr(sigma_expr) * I3
vm_expr = sqrt(1.5 * inner(sigma_dev, sigma_dev))

# ==========================================================
# 7. Loading Phase (100 Steps Ramp-Up)
# ==========================================================
load_steps = 100
print("\n" + "=" * 50)
print(">>> STARTING LOADING PHASE (100 Increments)")
print("=" * 50)

xdmf_load = XDMFFile(f"{file_name}_loading_phase.xdmf")
xdmf_load.parameters["flush_output"] = True
xdmf_load.parameters["functions_share_mesh"] = True

for step in range(1, load_steps + 1):
    factor = step / float(load_steps)
    load_factor.assign(factor)
    disp_factor.assign(factor) # Update the Constant value

    solve(F_res == 0, u, bcs=bcs, J=J_res, solver_parameters=solver_params)

    if step % 20 == 0 or step == load_steps:
        current_iop_mmHg = factor * iop_target_MPa * 7500.62
        print(f"Loading Step [{step:03d}/{load_steps}] -> Applied IOP = {current_iop_mmHg:.1f} mmHg")

        vm_func.assign(project(vm_expr, V_DG0))
        xdmf_load.write(u, float(step))
        xdmf_load.write(vm_func, float(step))

xdmf_load.close()
print("Loading phase completed successfully.\n")

# ==========================================================
# 8. Remodeling Phase (PT Subdomain 3 Only)
# ==========================================================
remodel_steps = 500
dt = 0.1
growth_rate = 0.005
theta_min = 0.5        # Max atrophy limit (50% residual volume)
gamma = 2.0
E1_crit = 0.02         # 2% tensile strain threshold

pt_tag = 3             # PT is Subdomain 3
pt_cells = np.where(sub_arr == pt_tag)[0]

print("=" * 50)
print(">>> STARTING REMODELING PHASE (PT Subdomain 3 Only)")
print("=" * 50)

xdmf_remodel = XDMFFile(f"{file_name}_remodeling_phase.xdmf")
xdmf_remodel.parameters["flush_output"] = True
xdmf_remodel.parameters["functions_share_mesh"] = True

theta_arr = np.array(theta_func.vector().get_local())
E1_arr = np.zeros(mesh.num_cells())

# Maintain full load throughout remodeling
load_factor.assign(1.0)
rad_disp_expr.factor = 1.0

for rstep in range(1, remodel_steps + 1):
    # 1. Project Green-Lagrange Strain Tensor
    E_tensor_proj = project(E_gl_expr, V_tensor)
    E_tensor_vals = E_tensor_proj.vector().get_local().reshape((-1, 3, 3))

    # 2. Extract Max Principal Strain (E1) via Eigenvalues in PT cells
    for cell_idx in pt_cells:
        dof = dof_map.cell_dofs(cell_idx)[0]
        E_mat = E_tensor_vals[cell_idx]
        
        # Symmetrize numerical tensor
        E_sym = 0.5 * (E_mat + E_mat.T)
        eigvals = np.linalg.eigvalsh(E_sym)
        E1_arr[cell_idx] = np.max(eigvals)

        # 3. Evolution ODE on PT cells
        E1_val = E1_arr[cell_idx]
        if E1_val > E1_crit:
            overload = (E1_val - E1_crit) / E1_crit
            curr_theta = theta_arr[dof]
            
            k_theta = growth_rate * ((curr_theta - theta_min) / (1.0 - theta_min)) ** gamma
            k_theta = max(k_theta, 0.0)

            d_theta = -k_theta * overload * dt * curr_theta
            theta_arr[dof] = np.clip(curr_theta + d_theta, theta_min, 1.0)

    # 4. Update Theta Field
    theta_func.vector().set_local(theta_arr)
    theta_func.vector().apply("insert")

    # 5. Re-solve Equilibrium with Updated Growth State
    solve(F_res == 0, u, bcs=bcs, J=J_res, solver_parameters=solver_params)

    # 6. Output & Logging
    if rstep % 5 == 0 or rstep == remodel_steps:
        pt_dofs = [dof_map.cell_dofs(c)[0] for c in pt_cells]
        avg_theta_pt = np.mean(theta_arr[pt_dofs])
        min_theta_pt = np.min(theta_arr[pt_dofs])
        max_E1_pt = np.max(E1_arr[pt_cells])

        print(f"[Remodeling Step {rstep:03d}/{remodel_steps}] "
              f"PT Avg Theta: {avg_theta_pt:.4f} | Min Theta: {min_theta_pt:.4f} | Max E1: {max_E1_pt:.4f}")

        E1_func.vector().set_local(E1_arr)
        E1_func.vector().apply("insert")
        vm_func.assign(project(vm_expr, V_DG0))

        xdmf_remodel.write(u, float(rstep))
        xdmf_remodel.write(theta_func, float(rstep))
        xdmf_remodel.write(E1_func, float(rstep))
        xdmf_remodel.write(vm_func, float(rstep))

xdmf_remodel.close()
print("\nRemodeling process finished! Results saved to XDMF.")