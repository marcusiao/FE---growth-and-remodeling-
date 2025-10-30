// Gmsh project created on Tue Oct 14 19:42:08 2025
SetFactory("OpenCASCADE");

SetFactory("OpenCASCADE");

// --------------------------------------------------
// 1. Define cube size
// --------------------------------------------------
lc = 1;           // mesh size (can adjust to refine)
L = 5;            // length in mm

// --------------------------------------------------
// 2. Create cube geometry
// --------------------------------------------------
Box(1) = {0, 0, 0, L, L, L};

// --------------------------------------------------
// 3. Define physical groups for boundaries and volume
// --------------------------------------------------
Physical Volume("Cube_Volume") = {1};


// --------------------------------------------------
// 4. Mesh settings
// --------------------------------------------------
Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;
Mesh 3;
//+
Physical Surface("clamp", 13) = {3};
//+
Physical Surface("load", 14) = {4};
