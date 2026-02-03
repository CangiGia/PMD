function inJoints
    include_global

J1 = Joint_struct;
J1.type = 'rev-rev'; % Q-A
J1.iPindex = 1;
J1.jPindex = 5;
J1.L = 0.45;

J2 = Joint_struct;
J2.type = 'rev-tran'; % O-B
J2.iPindex = 2;
J2.jPindex = 4;
J2.iUindex = 1;

Joints = [J1; J2];