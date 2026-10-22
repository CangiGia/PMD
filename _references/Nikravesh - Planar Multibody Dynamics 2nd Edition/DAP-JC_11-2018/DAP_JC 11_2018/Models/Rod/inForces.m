function inForces
    include_global

S1 = Force_struct;
S1.type = 'weight';  % include the weight

S2 = Force_struct;
S2.type = 'user'; 

Forces = [S1; S2];
