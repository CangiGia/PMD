function inJoints
    include_global

J1 = Joint_struct;
J1.type = 'tran'; % torso slider
J1.iPindex = 2;
J1.jPindex = 1;
J1.iUindex = 2;
J1.jUindex = 1;

J2 = Joint_struct;
J2.type = 'rev'; % neck pin 0
J2.iPindex = 3;
J2.jPindex = 4;

J3 = Joint_struct;
J3.type = 'rev'; % neck pin 1
J3.iPindex = 5;
J3.jPindex = 6;

J4 = Joint_struct;
J4.type = 'rev'; % neck pin 2
J4.iPindex = 7;
J4.jPindex = 8;

J5 = Joint_struct;
J5.type = 'rev'; % neck pin 3
J5.iPindex = 9;
J5.jPindex = 10;

J6 = Joint_struct;
J6.type = 'rev'; % neck pin 4
J6.iPindex = 11;
J6.jPindex = 12;

J7 = Joint_struct;
J7.type = 'rev'; % neck pin 5
J7.iPindex = 13;
J7.jPindex = 14;

J8 = Joint_struct;
J8.type = 'rev'; % neck pin 6
J8.iPindex = 15;
J8.jPindex = 16;

J9 = Joint_struct;
J9.type = 'rev'; % neck pin 7
J9.iPindex = 17;
J9.jPindex = 18;

Joints = [J1; J2; J3; J4; J5; J6; J7; J8; J9];
