! +++++ REV JOINT :: GROUND - Link-1
! ------------------------------------------------------- !
! ------------------------------------------------------- !
marker create marker=.double_compound_pendulum.ground.ground_link_1_reference_frame_ori &
    location=0.0, 0.0, 0.0 &
    orientation=0.0, 0.0, 0.0

marker create marker=.double_compound_pendulum.link_1.link_1_ground_reference_frame_ori &
    location=(LOC_RELATIVE_TO({0.0,0.0,0.0},.double_compound_pendulum.ground.ground_link_1_reference_frame_ori)) &
    orientation=(ORI_RELATIVE_TO({0.0,0.0,0.0},.double_compound_pendulum.ground.ground_link_1_reference_frame_ori))

constraint create joint revolute  &
	joint_name = .double_compound_pendulum.ground_link_1_joint  &
	friction_enabled = no  &
	i_marker_name = .double_compound_pendulum.ground.ground_link_1_reference_frame_ori &
	j_marker_name = .double_compound_pendulum.link_1.link_1_ground_reference_frame_ori
! ------------------------------------------------------- !
! ------------------------------------------------------- !

! +++++ REV JOINT :: Link-1 - Link-2
! ------------------------------------------------------- ! 
! ------------------------------------------------------- !
marker create marker=.double_compound_pendulum.link_1.link_1_link_2_reference_frame_ori &
    location=(LOC_RELATIVE_TO({0.0,0.0,0.0},.double_compound_pendulum.link_1.link_1_link_2_reference_frame)) &
    orientation=0.0, 0.0, 0.0

marker create marker=.double_compound_pendulum.link_2.link_2_link_1_reference_frame_ori &
    location=(LOC_RELATIVE_TO({0.0,0.0,0.0},.double_compound_pendulum.link_1.link_1_link_2_reference_frame_ori)) &
    orientation=(ORI_RELATIVE_TO({0.0,0.0,0.0},.double_compound_pendulum.link_1.link_1_link_2_reference_frame_ori))   

constraint create joint revolute  &
    joint_name = .double_compound_pendulum.link_1_link_2_joint  &
    friction_enabled = no  &
    i_marker_name = .double_compound_pendulum.link_1.link_1_link_2_reference_frame_ori &
    j_marker_name = .double_compound_pendulum.link_2.link_2_link_1_reference_frame_ori
! ------------------------------------------------------- !
! ------------------------------------------------------- !