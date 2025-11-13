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
load_steps = 10
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
# 5. Remodeling Phase (maintain load)
# ======================================================
growth_steps = 10



dx = Measure("dx", domain=mesh, metadata={"quadrature_degree": 2})

Fg_space = TensorFunctionSpace(mesh, "CG", 1)
Fg = Function(Fg_space)
Fg.assign(project(Identity(d), Fg_space))

# Keep load at final value
load_factor.assign(target_load)
pressure_vec = load_factor * pressure_base

xdmf_growth = XDMFFile(mesh.mpi_comm(), "u_growth_phase.xdmf")
xdmf_growth.parameters["flush_output"] = True
xdmf_growth.parameters["functions_share_mesh"] = True

xdmf_Fg = XDMFFile(mesh.mpi_comm(), "Fg_growth_phase.xdmf")
xdmf_Fg.parameters["flush_output"] = True
xdmf_Fg.parameters["functions_share_mesh"] = True

print("=== Growth and Remodeling Phase ===")

# --- 0. Setup material property fields per element ---
DG0 = FunctionSpace(mesh, "DG", 0)  # piecewise constant per cell
mu_field = Function(DG0)
lmbda_field = Function(DG0)

# Initial material property
mu_new = 1.3 * mu  # new value after remodeling
lmbda_new = 1.2 * lmbda

mu_field.assign(Constant(mu))
lmbda_field.assign(Constant(lmbda))

# Stress threshold for remodeling
stress_threshold = 1e-1

updated_cell_count = 0

# --- Remodeling + growth loop ---
for gstep in range(1, growth_steps + 1):

    # --- Compute total deformation ---
    F_current = I + grad(u)
    Fe = F_current * inv(Fg)

    # --- Compute stress per element ---
    # Approximate Cauchy stress (linear elastic for simplicity)
    sigma = 2 * mu_field * sym(Fe - I) + lmbda_field * tr(sym(Fe - I)) * I

# --- Project stress magnitude to DG0 to get one value per cell ---
    sigma_proj = project(sqrt(inner(sigma, sigma)), DG0)
    sigma_cell_vals = sigma_proj.vector().get_local()

    # --- Average stress over all cells ---
    avg_cell_stress = np.mean(sigma_cell_vals)
    print(f"[Step {step}] Avg cell stress = {avg_cell_stress:.6e}")

    # --- Compute stress magnitude per cell ---
    # Project stress magnitude onto DG0 for cell-wise evaluation
    sigma_mag = project(sqrt(inner(sigma, sigma)), DG0)

    # --- Remodeling: update material properties based on threshold ---
    sigma_array = sigma_mag.vector().get_local()
    mu_array = mu_field.vector().get_local()
    lmbda_array = lmbda_field.vector().get_local()

    for i in range(len(sigma_array)):
        if sigma_array[i] > stress_threshold:
            mu_array[i] = mu_new
            lmbda_array[i] = lmbda_new
            updated_cell_count += 1

    print("updated_cell_count : " , updated_cell_count)

    # Assign updated values back to Functions
    mu_field.vector().set_local(mu_array)
    mu_field.vector().apply("insert")
    lmbda_field.vector().set_local(lmbda_array)
    lmbda_field.vector().apply("insert")

   # --- Compute average material properties for monitoring ---
    avg_mu = np.mean(mu_field.vector().get_local())
    avg_lmbda = np.mean(lmbda_field.vector().get_local())
    print(f"[Step {step}] Avg mu = {avg_mu:.4f}, Avg lambda = {avg_lmbda:.4f}")

    # --- Solve updated variational problem with new material ---
    F_tot = F_current * inv(Fg)
    C = F_tot.T * F_tot
    Ic = tr(C)
    J = det(F_tot)
    psi = (mu_field/2)*(Ic - 3) - mu_field*ln(J) + (lmbda_field/2)*(ln(J))**2

    Pi_total = psi*dx - dot(pressure_vec, u)*ds_measure(pressure_id)
    F_res = derivative(Pi_total, u, v)
    J_res = derivative(F_res, u, du)
    solve(F_res == 0, u, bcs=[bc_clamp], J=J_res)

    # --- Save results ---
    xdmf_growth.write(u, float(gstep))
    xdmf_Fg.write(Fg, float(gstep))







# for gstep in range(1, growth_steps + 1):
#     # --- Update elastic deformation ---
#     F_current = I + grad(u)
#     Fe = F_current * inv(Fg)

#     # --- Make Fe differentiable for UFL derivative ---
#     Fe_var = variable(Fe)

#     # --- Define strain energy in terms of Fe_var ---
#     C_e = Fe_var.T * Fe_var
#     Ic = tr(C_e)
#     J = det(Fe_var)
#     psi = (mu/2)*(Ic - 3) - mu*ln(J) + (lmbda/2)*(ln(J))**2

#     # --- Compute Cauchy stress ---
#     sigma = (1/J) * Fe_var * diff(psi, Fe_var).T

#     # --- Extract stress components ---
#     sigma_xx = sigma[0, 0]
#     sigma_yy = sigma[1, 1]
#     sigma_zz = sigma[2, 2]

#     # --- Compute growth tensor increment ---
#     F_x = T_recip * (sigma_xx - sigma_x) / sigma_x + 1
#     F_y = T_recip * (sigma_yy - sigma_y) / sigma_y + 1
#     F_z = T_recip * (sigma_zz - sigma_z) / sigma_z + 1

#     #     # --- Stability cap: limit growth per step ---
#     # F_x = max(1.0, min(F_x, 1.1))   # clamp between 1.0 and 1.1
#     # F_y = max(1.0, min(F_y, 1.1))   # clamp between 1.0 and 1.1
#     # F_y = max(1.0, min(F_y, 1.1))   # clamp between 1.0 and 1.1

#     Fg_increment = as_tensor([[F_x, 0, 0],
#                             [0, F_y, 0],
#                             [0, 0, F_z]])

#     # --- Update total growth tensor ---
#     Fg.assign(project(Fg_increment * Fg, Fg_space))


#     # --- Total deformation after growth ---
#     F_tot = F_current * inv(Fg)
#     C = F_tot.T*F_tot
#     Ic = tr(C)
#     J = det(F_tot)
#     psi = (mu/2)*(Ic - 3) - mu*ln(J) + (lmbda/2)*(ln(J))**2

#     # --- Update variational problem ---
#     Pi_total_growth = psi*dx - dot(pressure_vec, u)*ds_measure(pressure_id)
#     F_res_growth = derivative(Pi_total_growth, u, v)
#     J_res_growth = derivative(F_res_growth, u, du)

#     # --- Solve ---
#     solve(F_res_growth == 0, u, bcs=[bc_clamp], J=J_res_growth)

#     # --- Save and monitor stress ---
#     if gstep % 10 == 0:
#         # Compute Cauchy stress tensor (approx.)
#         sigma = project(mu*(grad(u) + grad(u).T) + lmbda*tr(grad(u))*I, 
#                         TensorFunctionSpace(mesh, "CG", 1))
#         sigma_mag = sqrt(inner(sigma, sigma))
#         avg_stress = assemble(sigma_mag*dx) / assemble(1*dx)
#         print(f"[Growth Step {gstep}] Avg stress = {avg_stress:.6e}")


#     xdmf_growth.write(u, float(gstep))
#     if gstep % 20 == 0:
#         xdmf_Fg.write(Fg, float(gstep))

# xdmf_growth.close()
# xdmf_Fg.close()
# print("✅ Remodeling phase finished and saved.")

# ======================================================
# 6. Unloading Phase (gradual release from remodeled state)
# ======================================================
unload_steps = 10
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

