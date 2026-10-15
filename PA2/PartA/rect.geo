//+
SetFactory("OpenCASCADE");
//+
Point(1) = {-0, -0, 0, 0.03};
//+
Point(2) = {2, -0, 0, 0.03};
//+
Point(3) = {2, 1, 0, 0.03};
//+
Point(4) = {-0, 1, 0, 0.03};
//+
Line(1) = {1, 1};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {1, 4};
//+
Physical Curve("pec", 1000000) = {4, 1, 2, 3};
//+
Curve Loop(1) = {3, -4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Physical Surface("dielectric1", 1) = {1};
