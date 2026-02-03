function inJoints
    include_global

J1 = Joint_struct;
J1.type = 'tran';
J1.iPindex = 2;
J1.jPindex = 1;
J1.iUindex = 2;
J1.jUindex = 1;

J2 = Joint_struct;
J2.type = 'rev';
J2.iPindex = 2;
J2.jPindex = 3;

Joints = [J1; J2];
