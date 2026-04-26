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
zero = Constant((0.0, 0.0))
clamp_id = 13
bc_clamp = DirichletBC(V, zero, boundaries, clamp_id)

E = 1.0
nu = 0.3
mu = Constant(E / (2*(1+nu)))
lmbda = Constant(E*nu / ((1+nu)*(1-2*nu)))

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
target_load = 1e-1
pressure_id = 14
pressure_base = as_vector((0, -1))
load_factor = Constant(0.0)
pressure_vec = load_factor * pressure_base
ds_measure = Measure("ds", domain=mesh, subdomain_data=boundaries)

F = I + grad(u)
C = F.T * F
Ic = tr(C)
J = det(F)
psi = (mu / 2) * (Ic - 3) - mu * ln(J) + (lmbda / 2) * (ln(J))**2

Pi_total = psi * dx_c - dot(pressure_vec, u) * ds_measure(pressure_id)
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
    load_factor.assign(step / load_steps * target_load)
    solve(F_res == 0, u, bcs=[bc_clamp], J=J_res, solver_parameters=solver_params)
    xdmf_load.write(u, float(step))
xdmf_load.close()
print("✅ Loading phase finished and saved.")

# ======================================================
# 5. Remodeling Phase (maintain load)
# ======================================================
growth_steps = 100
dt_growth = 0.1
growth_rate = 0.05

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

# --- Define variational problem once outside the loop ---
F_current_sym = I + grad(u)
F_tot_sym = F_current_sym * inv(Fg)
C_sym = F_tot_sym.T * F_tot_sym
Ic_sym = tr(C_sym)
J_sym = det(F_tot_sym)
psi_growth = (mu/2)*(Ic_sym - 3) - mu*ln(J_sym) + (lmbda/2)*(ln(J_sym))**2
Pi_total_growth = psi_growth*dx_c - dot(pressure_vec, u)*ds_measure(pressure_id)
F_res_growth = derivative(Pi_total_growth, u, v)
J_res_growth = derivative(F_res_growth, u, du)

print("=== Remodeling Phase ===")
for gstep in range(1, growth_steps + 1):
    # --- Update elastic deformation using current u ---
    F_current = I + grad(u)
    Fe = F_current * inv(Fg)
    E_tensor = 0.5*(Fe.T*Fe - I)
    strain_mag = sqrt(inner(E_tensor, E_tensor))

    # --- Growth increment ---
    # Optimized: combine growth multiplier and previous Fg into one projection
    Fg.assign(project((1.0 + growth_rate*strain_mag*dt_growth) * Fg, Fg_space))

    # --- Solve ---
    solve(F_res_growth == 0, u, bcs=[bc_clamp], J=J_res_growth, solver_parameters=solver_params)

    # --- Save ---
    if gstep % 10 == 0:
        avg_strain = assemble(strain_mag*dx_c)/assemble(1*dx_c)
        print(f"[Growth Step {gstep}] Avg strain = {avg_strain:.6e}")

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
F_tot_un = F_current_un * inv(Fg)
C_un = F_tot_un.T * F_tot_un
Ic_un = tr(C_un)
J_un = det(F_tot_un)
psi_un = (mu / 2) * (Ic_un - 3) - mu * ln(J_un) + (lmbda / 2) * (ln(J_un))**2
Pi_total_unload = psi_un * dx_c - dot(pressure_vec, u) * ds_measure(pressure_id)
F_res_unload = derivative(Pi_total_unload, u, v)
J_res_unload = derivative(F_res_unload, u, du)

print("=== Unloading Phase ===")
for step in range(1, unload_steps + 1):
    # Gradually reduce load
    load_factor.assign(target_load * (1 - step / unload_steps))

    # --- Solve equilibrium ---
    solve(F_res_unload == 0, u, bcs=[bc_clamp], J=J_res_unload, solver_parameters=solver_params)

    # --- Save displacement ---
    xdmf_unload.write(u, float(step))

xdmf_unload.close()
print("✅ Unloading phase finished and saved.")
