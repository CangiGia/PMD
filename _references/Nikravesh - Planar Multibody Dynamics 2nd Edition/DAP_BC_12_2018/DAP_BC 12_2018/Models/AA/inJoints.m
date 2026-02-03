function inJoints
    include_global

J1 = Joint_struct;
J1.type = 'rev'; %Q
J1.iPindex = 1;
J1.jPindex = 8;

J2 = Joint_struct;
J2.type = 'rev'; % A
J2.iPindex = 2;
J2.jPindex = 3;

J3 = Joint_struct;
J3.type = 'rev'; % B
J3.iPindex = 4;
J3.jPindex = 5;

J4 = Joint_struct;
J4.type = 'rev'; % O
J4.iPindex = 6;
J4.jPindex = 7;

Joints = [J1; J2; J3; J4];
