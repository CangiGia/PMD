function inJoints
    include_global

J1 = Joint_struct;
J1.type = 'tran';
J1.iPindex = 2;
J1.jPindex = 1;
J1.iUindex = 2;
J1.jUindex = 1;

Joints = [J1];
