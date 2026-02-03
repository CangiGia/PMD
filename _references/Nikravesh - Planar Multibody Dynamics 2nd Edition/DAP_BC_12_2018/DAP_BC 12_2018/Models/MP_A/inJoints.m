function inJoints
    include_global

J1 = Joint_struct;
J1.type = 'rev'; %A
J1.iPindex = 1;
J1.jPindex = 7;

J2 = Joint_struct;
J2.type = 'rev'; % Q
J2.iPindex = 5;
J2.jPindex = 6;

J3 = Joint_struct;
J3.type = 'rev'; % O
J3.iPindex = 4;
J3.jPindex = 8;

J4 = Joint_struct;
J4.type = 'tran'; % O-B
J4.iUindex = 1;
J4.jUindex = 2;
J4.iPindex = 2;
J4.jPindex = 8;

Joints = [J1; J2; J3; J4];
