function inUvectors
    include_global

V1 = Unit_struct;
V1.Bindex = 2;
V1.ulocal  = [0.47; -0.88];

V2 = Unit_struct;
V2.Bindex = 3;
V2.ulocal  = [1; 0];

Uvectors = [V1; V2];
