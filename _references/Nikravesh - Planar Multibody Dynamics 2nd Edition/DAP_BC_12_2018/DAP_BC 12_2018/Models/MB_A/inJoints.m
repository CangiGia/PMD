function inJoints
    include_global

J1 = Joint_struct;
J1.type = 'rev'; % O
J1.iPindex = 1;
J1.jPindex = 2;

J2 = Joint_struct;
J2.type = 'rev'; % A
J2.iPindex = 3;
J2.jPindex = 4;

J3 = Joint_struct;
J3.type = 'rev-tran'; % front fork D
J3.iPindex = 7;
J3.jPindex = 8;
J3.iUindex = 1;
J3.L = 0;

Joints = [J1; J2; J3];
