function inForces
    include_global

S1 = Force_struct;
S1.type = 'rot_sda'; % torsional spring-damper at S1 (no gravity in orbit)
S1.iBindex = 2;
S1.jBindex = 3;
S1.k = 30;
S1.theta0 = -1.08219905;
S1.dc = 7.5;
S1.T_a = 0;

Forces = [S1];
