function inJoints
    include_global

J1 = Joint_struct;
J1.type = 'rev'; % A
J1.iPindex = 1;
J1.jPindex = 2;

J2 = Joint_struct;
J2.type = 'rev'; % G
J2.iPindex = 3;
J2.jPindex = 4;

J3 = Joint_struct;
J3.type = 'rev'; % D
J3.iPindex = 5;
J3.jPindex = 6;

J4 = Joint_struct;
J4.type = 'rev'; % C
J4.iPindex = 7;
J4.jPindex = 8;

J5 = Joint_struct;
J5.type = 'rev'; % I
J5.iPindex = 9;
J5.jPindex = 10;

J6 = Joint_struct;
J6.type = 'rev'; % E
J6.iPindex = 11;
J6.jPindex = 12;

J7 = Joint_struct;
J7.type = 'rev'; % K
J7.iPindex = 13;
J7.jPindex = 14;

Joints = [J1; J2; J3; J4; J5; J6; J7];
