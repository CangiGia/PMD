function inJoints
    include_global

J1 = Joint_struct;
J1.type = 'rev'; % A
J1.iPindex = 3;
J1.jPindex = 1;

J2 = Joint_struct;
J2.type = 'rev'; % C
J2.iPindex = 4;
J2.jPindex = 5;

J3 = Joint_struct;
J3.type = 'rev'; % D
J3.iPindex = 6;
J3.jPindex = 8;

J4 = Joint_struct;
J4.type = 'rev'; % B
J4.iPindex = 10;
J4.jPindex = 2;

J5 = Joint_struct;
J5.type = 'rev'; % E
J5.iPindex = 9;
J5.jPindex = 11;

J6 = Joint_struct;
J6.type = 'rev'; % F
J6.iPindex = 12;
J6.jPindex = 13;

J7 = Joint_struct;
J7.type = 'rev'; % K
J7.iPindex = 7;
J7.jPindex = 14;

J8 = Joint_struct;
J8.type = 'rel-tran'; % hydraulic actuator A-D
J8.iPindex = 6;
J8.jPindex = 1;
J8.iBindex = 2;
J8.jBindex = 0;
J8.iFunct = 1;

Joints = [J1; J2; J3; J4; J5; J6; J7; J8];
