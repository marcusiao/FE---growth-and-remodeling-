import meshio
from dolfin import *
import numpy as np
import os

# ======================================================
# 1. Setup 
# ======================================================
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# ======================================================
# 2. Load mesh and define spaces
# ======================================================
mesh = Mesh("Job1.xml")
facet_markers = MeshFunction("size_t", mesh, "Job1_facet_region.xml")

V = VectorFunctionSpace(mesh, "CG", 1)

# ======================================================
# 3. Boundary conditions and material parameters
# ======================================================
zero = Constant((0.0, 0.0, 0.0))
clamp_id = 1
bc_clamp = DirichletBC(V, zero, facet_markers, clamp_id)

E = 1.0
nu = 0.3
mu = Constant(E / (2*(1+nu)))
lmbda = Constant(E*nu / ((1+nu)*(1-2*nu)))

u = Function(V)
du = TrialFunction(V)
v = TestFunction(V)

d = mesh.geometry().dim()
I = Identity(d)

# ======================================================
# 4. Loading Phase (step load)
# ======================================================
load_steps = 100
target_load = 5e-2
pressure_id = 2
pressure_base = as_vector((0, -1, 0.0))
load_factor = Constant(0.0)
pressure_vec = load_factor * pressure_base
ds_measure = Measure("ds", domain=mesh, subdomain_data=facet_markers)

F = I + grad(u)
C = F.T*F
Ic = tr(C)
J = det(F)
psi = (mu/2)*(Ic - 3) - mu*ln(J) + (lmbda/2)*(ln(J))**2

Pi_total = psi*dx - dot(pressure_vec, u)*ds_measure(pressure_id)
F_res = derivative(Pi_total, u, v)
J_res = derivative(F_res, u, du)

xdmf_load = XDMFFile(mesh.mpi_comm(), "u_loading_phase.xdmf")
xdmf_load.parameters["flush_output"] = True
xdmf_load.parameters["functions_share_mesh"] = True

print("=== Loading Phase ===")
for step in range(1, load_steps + 1):
    load_factor.assign(step / load_steps * target_load)
    solve(F_res == 0, u, bcs=[bc_clamp], J=J_res)
    xdmf_load.write(u, float(step))
xdmf_load.close()
print("✅ Loading phase finished and saved.")

# ======================================================
# ======================================================
# ======================================================
# 5. Growth + Remodeling Phase (memory-safe version)
# ======================================================
growth_steps = 100

# --- LT2 Growth Law Parameters ---
sigma_x = 3.0   # threshold stress (x)
sigma_y = 3.0   # threshold stress (y)
sigma_z = 3.0   # threshold stress (z)
T_recip = 0.001 # time scaling for growth rate

dx = Measure("dx", domain=mesh, metadata={"quadrature_degree": 2})
print("⚙️ Using quadrature_degree = 2 to avoid large memory allocations")

# --- Growth tensor setup ---
Fg_space = TensorFunctionSpace(mesh, "CG", 1)
Fg = Function(Fg_space)
Fg.assign(project(Identity(d), Fg_space))

# --- Keep load constant at final value ---
load_factor.assign(target_load)
pressure_vec = load_factor * pressure_base

xdmf_growth = XDMFFile(mesh.mpi_comm(), "u_growth_phase.xdmf")
xdmf_growth.parameters["flush_output"] = True
xdmf_growth.parameters["functions_share_mesh"] = True

xdmf_Fg = XDMFFile(mesh.mpi_comm(), "Fg_growth_phase.xdmf")
xdmf_Fg.parameters["flush_output"] = True
xdmf_Fg.parameters["functions_share_mesh"] = True

print("=== Growth and Remodeling Phase ===")

# --- Material properties as piecewise constants ---
DG0 = FunctionSpace(mesh, "DG", 0)
mu_field = Function(DG0)
lmbda_field = Function(DG0)

mu_field.assign(Constant(mu))
lmbda_field.assign(Constant(lmbda))

#   define the trigger threshold and the updated material properties 
mu_new = 1.3 * mu
lmbda_new = 1.2 * lmbda
stress_threshold = 1e-1

for gstep in range(1, growth_steps + 1):

    #  monitor how many cells are being updated at each step
    updated_cell_count = 0

    # --- Elastic deformation ---
    F_current = I + grad(u)
    Fe = F_current * inv(Fg)

    # --- Project Fe to a low-order space for per-cell stress computation ---
    Vdg = TensorFunctionSpace(mesh, "DG", 0)
    Fe_proj = project(Fe, Vdg, form_compiler_parameters={"quadrature_degree": 2})

    # --- Approximate Cauchy stress (no symbolic differentiation) to cut down computational cost---
    # σ ≈ μ (Fe_proj + Fe_projᵀ − 2I) + λ tr(Fe_proj − I) I
    sigma_approx = project(
        mu_field * (Fe_proj + Fe_proj.T - 2 * I) +
        lmbda_field * tr(Fe_proj - I) * I,
        Vdg, form_compiler_parameters={"quadrature_degree": 2}
    )

    # --- Stress magnitude per cell ---
    sigma_mag = project(sqrt(inner(sigma_approx, sigma_approx)), DG0,
                        form_compiler_parameters={"quadrature_degree": 2})
    sigma_array = sigma_mag.vector().get_local()

    # --- Remodeling rule: update material where stress exceeds threshold ---
    mu_array = mu_field.vector().get_local()
    lmbda_array = lmbda_field.vector().get_local()

    for i in range(len(sigma_array)):
        if sigma_array[i] > stress_threshold:
            mu_array[i] = mu_new
            lmbda_array[i] = lmbda_new
            updated_cell_count += 1

#   load the new value into the model
    mu_field.vector().set_local(mu_array)
    mu_field.vector().apply("insert")
    lmbda_field.vector().set_local(lmbda_array)
    lmbda_field.vector().apply("insert")

#  calculate the average material properties for monitoring
    avg_mu = np.mean(mu_field.vector().get_local())
    avg_lmbda = np.mean(lmbda_field.vector().get_local())
    avg_sigma = np.mean(sigma_array)

    print(f"[Step {gstep}] Avg σ = {avg_sigma:.4e}, Avg μ = {avg_mu:.4f}, "
          f"Avg λ = {avg_lmbda:.4f}, Updated cells = {updated_cell_count}")

    # --- Growth update ---
    sigma_xx = sigma_approx[0, 0]
    sigma_yy = sigma_approx[1, 1]
    sigma_zz = sigma_approx[2, 2]


####  the LT2 growth law equations
    F_x = T_recip * (sigma_xx - sigma_x) / sigma_x + 1
    F_y = T_recip * (sigma_yy - sigma_y) / sigma_y + 1
    F_z = T_recip * (sigma_zz - sigma_z) / sigma_z + 1

    Fg_increment = as_tensor([[F_x, 0, 0],
                              [0, F_y, 0],
                              [0, 0, F_z]])

    Fg.assign(project(Fg_increment * Fg, Fg_space,
                      form_compiler_parameters={"quadrature_degree": 2}))

    # --- Solve updated equilibrium ---
    F_tot = F_current * inv(Fg)
    C = F_tot.T * F_tot
    Ic = tr(C)
    J = det(F_tot)
    psi = (mu_field / 2) * (Ic - 3) - mu_field * ln(J) + (lmbda_field / 2) * (ln(J)) ** 2

    Pi_total_growth = psi * dx - dot(pressure_vec, u) * ds_measure(pressure_id)
    F_res_growth = derivative(Pi_total_growth, u, v)
    J_res_growth = derivative(F_res_growth, u, du)

    solve(F_res_growth == 0, u, bcs=[bc_clamp], J=J_res_growth)

    # --- Output ---
    xdmf_growth.write(u, float(gstep))
    if gstep % 20 == 0:
        xdmf_Fg.write(Fg, float(gstep))

xdmf_growth.close()
xdmf_Fg.close()
print("✅ Growth + Remodeling phase finished and saved.")


# ======================================================
# 6. Unloading Phase (gradual release from remodeled state)
# ======================================================
unload_steps = 100
xdmf_unload = XDMFFile(mesh.mpi_comm(), "u_unloading_phase.xdmf")
xdmf_unload.parameters["flush_output"] = True
xdmf_unload.parameters["functions_share_mesh"] = True

print("=== Unloading Phase ===")
for step in range(1, unload_steps + 1):
    # Gradually reduce load
    load_factor.assign(target_load * (1 - step / unload_steps))
    pressure_vec = load_factor * pressure_base

    # --- Compute total deformation including remodeling ---
    F_current = I + grad(u)
    F_tot = F_current * inv(Fg)  # important: include growth tensor
    C = F_tot.T * F_tot
    Ic = tr(C)
    J = det(F_tot)
    psi = (mu / 2) * (Ic - 3) - mu * ln(J) + (lmbda / 2) * (ln(J))**2

    # --- Define variational problem with updated load ---
    Pi_total_unload = psi * dx - dot(pressure_vec, u) * ds_measure(pressure_id)
    F_res_unload = derivative(Pi_total_unload, u, v)
    J_res_unload = derivative(F_res_unload, u, du)

    # --- Solve equilibrium ---
    solve(F_res_unload == 0, u, bcs=[bc_clamp], J=J_res_unload)

    # --- Save displacement ---
    xdmf_unload.write(u, float(step))

xdmf_unload.close()
print("✅ Unloading phase finished and saved.")

