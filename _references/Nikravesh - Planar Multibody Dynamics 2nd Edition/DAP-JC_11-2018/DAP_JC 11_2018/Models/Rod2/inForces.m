function inForces
    include_global

F1 = Force_struct;
F1.type = 'user'; 

F2 = Force_struct;
F2.type = 'weight';  % include the weight

Forces = [F1; F2];