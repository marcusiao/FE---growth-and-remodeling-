import os
import numpy as np
import meshio

os.chdir(os.path.dirname(os.path.realpath(__file__)))

file_name = "cube"
inp = meshio.read(f"{file_name}.inp")

# ==========================================
# 1. Generate Domain Pair with Material Sets
# ==========================================
tet_cells = inp.get_cells_type("tetra")
tet_data = np.zeros(len(tet_cells), dtype=np.int32)

print("Found Abaqus Cell Sets:", list(inp.cell_sets.keys()))

# Assign Subdomain IDs based on Volume Sets:
# ID = 1 -> SOFT, ID = 2 -> HARD
for set_name, cell_blocks in inp.cell_sets.items():
    set_upper = set_name.upper()
    mat_id = None
    
    if "SOFT" in set_upper:
        mat_id = 1
    elif "HARD" in set_upper:
        mat_id = 2
        
    if mat_id is not None:
        for block in cell_blocks:
            for elem_idx in np.atleast_1d(block):
                elem_idx = int(elem_idx)
                if elem_idx < len(tet_cells):
                    tet_data[elem_idx] = mat_id

meshio.write_points_cells(
    f"{file_name}_domain.xdmf",
    inp.points,
    [("tetra", tet_cells)],
    cell_data={"name_to_read": [tet_data]}
)
print(f"Generated {file_name}_domain.xdmf and .h5 with subdomains 1 (soft) and 2 (hard)!")


# ==========================================
# 2. Extract Boundary Triangles & Map Surfaces
# ==========================================
TET_FACE_MAP = {
    "S1": [0, 1, 2],
    "S2": [0, 3, 1],
    "S3": [1, 3, 2],
    "S4": [2, 3, 0],
}

triangles = []
tri_data = []

for set_name, cell_blocks in inp.cell_sets.items():
    set_upper = set_name.upper()
    
    marker_id = None
    if "FIXED" in set_upper:
        marker_id = 1
    elif "LOAD" in set_upper:
        marker_id = 2
        
    if marker_id is not None:
        face_id = None
        for s_key in TET_FACE_MAP:
            if set_upper.endswith(s_key):
                face_id = s_key
                break
                
        if face_id is not None:
            local_nodes = TET_FACE_MAP[face_id]
            for block in cell_blocks:
                for elem_idx in np.atleast_1d(block):
                    elem_idx = int(elem_idx)
                    if elem_idx < len(tet_cells):
                        tet = tet_cells[elem_idx]
                        tri_face = [tet[local_nodes[0]], tet[local_nodes[1]], tet[local_nodes[2]]]
                        triangles.append(tri_face)
                        tri_data.append(marker_id)


# ==========================================
# 3. Generate Facet Region Pair (2D Triangles)
# ==========================================
if len(triangles) > 0:
    tri_cells = np.array(triangles, dtype=np.int64)
    tri_data_arr = np.array(tri_data, dtype=np.int32)

    meshio.write_points_cells(
        f"{file_name}_facet_region.xdmf",
        inp.points,
        [("triangle", tri_cells)],
        cell_data={"name_to_read": [tri_data_arr]}
    )
    print(f"Generated {file_name}_facet_region.xdmf and .h5 successfully!")
else:
    print("Warning: No matching surface faces found in cell sets!")