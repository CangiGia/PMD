function inUvectors
    include_global

V1 = Unit_struct;
V1.Bindex = 0;
V1.ulocal  = [1; 0];

V2 = Unit_struct;
V2.Bindex = 1;
V2.ulocal  = [1; 0];

Uvectors = [V1; V2];
