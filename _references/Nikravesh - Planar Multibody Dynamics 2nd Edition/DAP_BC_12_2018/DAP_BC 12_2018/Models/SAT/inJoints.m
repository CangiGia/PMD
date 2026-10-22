function inJoints
    include_global

J1 = Joint_struct;
J1.type = 'rev'; % B
J1.iPindex = 4;
J1.jPindex = 1;

J2 = Joint_struct;
J2.type = 'rev'; % A
J2.iPindex = 5;
J2.jPindex = 6;

J3 = Joint_struct;
J3.type = 'rev'; % S1 (links A-S1, S1-C)
J3.iPindex = 7;
J3.jPindex = 8;

J4 = Joint_struct;
J4.type = 'rev'; % S1 (links S1-C, S1-L1)
J4.iPindex = 9;
J4.jPindex = 11;

J5 = Joint_struct;
J5.type = 'rev'; % L1 (unlocked)
J5.iPindex = 12;
J5.jPindex = 13;

J6 = Joint_struct;
J6.type = 'rev'; % C
J6.iPindex = 10;
J6.jPindex = 2;

J7 = Joint_struct;
J7.type = 'rev'; % D
J7.iPindex = 14;
J7.jPindex = 3;

Joints = [J1; J2; J3; J4; J5; J6; J7];
