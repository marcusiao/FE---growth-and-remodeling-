import os
import numpy as np
import meshio

os.chdir(os.path.dirname(os.path.realpath(__file__)))

file_name = "3D_ONH_corase"
inp = meshio.read(f"{file_name}.inp")

# ==========================================================
# 1. Generate Domain Pair with Material Subdomains (Volume)
# ==========================================================
tet_cells = inp.get_cells_type("tetra")
tet_data = np.zeros(len(tet_cells), dtype=np.int32)

print("Found Abaqus Cell Sets:", list(inp.cell_sets.keys()))

# Subdomain Material Mapping:
# 1: CHOROID | 2: LC | 3: PT | 4: SCLERA
for set_name, cell_blocks in inp.cell_sets.items():
    set_upper = set_name.upper()
    mat_id = None
    
    if "CHOROID" in set_upper:
        mat_id = 1
    elif "LC" in set_upper:
        mat_id = 2
    elif "PT" in set_upper:
        mat_id = 3
    elif "SCLERA" in set_upper:
        mat_id = 4
        
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
print(f"Generated {file_name}_domain.xdmf and .h5 with material IDs (1: Choroid, 2: LC, 3: PT, 4: Sclera)!")


# ==========================================================
# 2. Extract Boundary Triangles & Map Surfaces (Facets)
# ==========================================================
# Abaqus C3D4 local face definitions (0-indexed node mappings)
# S1=(1,2,3), S2=(1,4,2), S3=(2,4,3), S4=(3,4,1)
TET_FACE_MAP = {
    "S1": [0, 1, 2],
    "S2": [0, 3, 1],
    "S3": [1, 3, 2],
    "S4": [2, 3, 0],
}

triangles = []
tri_data = []

# Boundary Facet Mapping:
# 1: SYMMETRY | 2: IOP | 3: FIXED | 4: OUT
for set_name, cell_blocks in inp.cell_sets.items():
    set_upper = set_name.upper()
    
    marker_id = None
    if "SYMMETRY" in set_upper:
        marker_id = 1
    elif "IOP" in set_upper:
        marker_id = 2
    elif "FIXED" in set_upper:
        marker_id = 3
    elif "OUT" in set_upper:
        marker_id = 4
        
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


# ==========================================================
# 3. Generate Facet Region Pair (2D Triangles)
# ==========================================================
if len(triangles) > 0:
    tri_cells = np.array(triangles, dtype=np.int64)
    tri_data_arr = np.array(tri_data, dtype=np.int32)

    meshio.write_points_cells(
        f"{file_name}_facet_region.xdmf",
        inp.points,
        [("triangle", tri_cells)],
        cell_data={"name_to_read": [tri_data_arr]}
    )
    print(f"Generated {file_name}_facet_region.xdmf and .h5 successfully with {len(tri_cells)} boundary facets!")
else:
    print("Warning: No matching surface faces found in cell sets!")