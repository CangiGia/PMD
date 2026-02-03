function inJoints
    include_global

J1 = Joint_struct;
J1.type = 'rev';
J1.iPindex = 1;
J1.jPindex = 3;

J2 = Joint_struct;
J2.type = 'rev';
J2.iPindex = 2;
J2.jPindex = 4;

J3 = Joint_struct;
J3.type = 'disc';
J3.iBindex = 2;
J3.R = 0.1;
J3.x0 = 0.2;

J4 = Joint_struct;
J4.type = 'disc';
J4.iBindex = 3;
J4.R = 0.1;
J4.x0 = 0.8;
J4.p0 = 0; % default

J5 = Joint_struct;
J5.type = 'rel-rot'; % motor driver
J5.iBindex = 2;
J5.jBindex = 1;
J5.iFunct = 1;

Joints = [J1; J2; J3; J4; J5];
