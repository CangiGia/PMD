function inForces
    include_global

F1 = Force_struct;
F1.type = 'weight'; 

F2 = Force_struct;
F2.type = 'user'; 

Forces = [F1; F2];
