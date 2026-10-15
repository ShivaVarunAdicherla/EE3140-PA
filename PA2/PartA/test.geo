// Gmsh project created on Fri Oct  9 18:06:01 2026
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 5, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Plane Surface(1) = {1};
//+
Physical Curve("pec", 1000000) = {1};
//+
Physical Surface("dielectric1", 1) = {1};
