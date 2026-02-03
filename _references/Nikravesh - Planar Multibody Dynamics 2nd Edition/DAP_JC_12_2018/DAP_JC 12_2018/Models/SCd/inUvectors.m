function inUvectors
    include_global

U1 = Unit_struct;
U1.Bindex = 0;
U1.ulocal  = [ 1.0; 0];

U2 = Unit_struct;
U2.Bindex = 3;
U2.ulocal  = [ 1.0; 0];

Uvectors = [U1; U2];
