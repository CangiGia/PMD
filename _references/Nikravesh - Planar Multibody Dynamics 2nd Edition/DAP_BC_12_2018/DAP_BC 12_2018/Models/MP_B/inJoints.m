function inJoints
    include_global

J1 = Joint_struct;
J1.type = 'rev'; %A
J1.iPindex = 1;
J1.jPindex = 7;

J2 = Joint_struct;
J2.type = 'rev-tran'; % O-B
J2.iPindex = 2;
J2.jPindex = 4;
J2.iUindex = 1;

J3 = Joint_struct;
J3.type = 'rev'; % Q
J3.iPindex = 5;
J3.jPindex = 6;

Joints = [J1; J2; J3];
