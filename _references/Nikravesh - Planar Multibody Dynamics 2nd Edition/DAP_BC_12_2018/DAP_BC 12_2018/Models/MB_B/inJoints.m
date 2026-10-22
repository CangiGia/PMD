function inJoints
    include_global

J1 = Joint_struct;
J1.type = 'rev'; % O
J1.iPindex = 1;
J1.jPindex = 2;

J2 = Joint_struct;
J2.type = 'rev'; % A (chainstay-seatstay)
J2.iPindex = 3;
J2.jPindex = 4;

J3 = Joint_struct;
J3.type = 'rev'; % B
J3.iPindex = 7;
J3.jPindex = 8;

J4 = Joint_struct;
J4.type = 'rev'; % Q
J4.iPindex = 9;
J4.jPindex = 10;

J5 = Joint_struct;
J5.type = 'rev'; % A (rear wheel)
J5.iPindex = 5;
J5.jPindex = 6;

J6 = Joint_struct;
J6.type = 'rev-tran'; % front fork E
J6.iPindex = 13;
J6.jPindex = 14;
J6.iUindex = 1;
J6.L = 0;

Joints = [J1; J2; J3; J4; J5; J6];
