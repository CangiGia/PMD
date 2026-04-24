! +++++ REV JOINT :: GROUND - Link-1
! ------------------------------------------------------- !
! ------------------------------------------------------- !
marker create marker=.four_bar_linkage.ground.ground_link_1_reference_frame_ori &
    location=0.0, 0.0, 0.0 &
    orientation=0.0, 0.0, 0.0

marker create marker=.four_bar_linkage.link_1.link_1_ground_reference_frame_ori &
    location=(LOC_RELATIVE_TO({0.0,0.0,0.0},.four_bar_linkage.ground.ground_link_1_reference_frame_ori)) &
    orientation=(ORI_RELATIVE_TO({0.0,0.0,0.0},.four_bar_linkage.ground.ground_link_1_reference_frame_ori))

constraint create joint revolute  &
	joint_name = .four_bar_linkage.ground_link_1_joint  &
	friction_enabled = no  &
	i_marker_name = .four_bar_linkage.ground.ground_link_1_reference_frame_ori &
	j_marker_name = .four_bar_linkage.link_1.link_1_ground_reference_frame_ori
! ------------------------------------------------------- !
! ------------------------------------------------------- !

! +++++ REV JOINT :: Link-1 - Link-2
! ------------------------------------------------------- ! 
! ------------------------------------------------------- !
marker create marker=.four_bar_linkage.link_1.link_1_link_2_reference_frame_ori &
    location=(LOC_RELATIVE_TO({0.0,0.0,0.0},.four_bar_linkage.link_1.link_1_link_2_reference_frame)) &
    orientation=0.0, 0.0, 0.0

marker create marker=.four_bar_linkage.link_2.link_2_link_1_reference_frame_ori &
    location=(LOC_RELATIVE_TO({0.0,0.0,0.0},.four_bar_linkage.link_1.link_1_link_2_reference_frame_ori)) &
    orientation=(ORI_RELATIVE_TO({0.0,0.0,0.0},.four_bar_linkage.link_1.link_1_link_2_reference_frame_ori))   

constraint create joint revolute  &
    joint_name = .four_bar_linkage.link_1_link_2_joint  &
    friction_enabled = no  &
    i_marker_name = .four_bar_linkage.link_1.link_1_link_2_reference_frame_ori &
    j_marker_name = .four_bar_linkage.link_2.link_2_link_1_reference_frame_ori
! ------------------------------------------------------- !
! ------------------------------------------------------- !

! +++++ REV JOINT :: Link-2 - Link-3
! ------------------------------------------------------- ! 
! ------------------------------------------------------- !
marker create marker=.four_bar_linkage.link_2.link_2_link_3_reference_frame_ori &
    location=(LOC_RELATIVE_TO({0.0,0.0,0.0},.four_bar_linkage.link_2.link_2_link_3_reference_frame)) &
    orientation=0.0, 0.0, 0.0

marker create marker=.four_bar_linkage.link_3.link_3_link_2_reference_frame_ori &
    location=(LOC_RELATIVE_TO({0.0,0.0,0.0},.four_bar_linkage.link_2.link_2_link_3_reference_frame_ori)) &
    orientation=(ORI_RELATIVE_TO({0.0,0.0,0.0},.four_bar_linkage.link_2.link_2_link_3_reference_frame_ori))   

constraint create joint revolute  &
    joint_name = .four_bar_linkage.link_2_link_3_joint  &
    friction_enabled = no  &
    i_marker_name = .four_bar_linkage.link_2.link_2_link_3_reference_frame_ori &
    j_marker_name = .four_bar_linkage.link_3.link_3_link_2_reference_frame_ori
! ------------------------------------------------------- !
! ------------------------------------------------------- !

! +++++ REVOLUTE JOINT :: Link-3 - GROUND
! ------------------------------------------------------- !
! ------------------------------------------------------- !
marker create marker=.four_bar_linkage.link_3.link_3_ground_reference_frame_ori &
    location=(LOC_RELATIVE_TO({0.0,0.0,0.0},.four_bar_linkage.link_3.link_3_ground_reference_frame)) &
    orientation=0.0, 0.0, 0.0

marker create marker=.four_bar_linkage.ground.ground_link_3_reference_frame_ori &
    location=(LOC_RELATIVE_TO({0.0,0.0,0.0},.four_bar_linkage.link_3.link_3_ground_reference_frame_ori)) &
    orientation=(ORI_RELATIVE_TO({0.0,0.0,0.0},.four_bar_linkage.link_3.link_3_ground_reference_frame_ori))

constraint create joint revolute &
    joint_name = .four_bar_linkage.link_3_ground_joint &
    i_marker_name = .four_bar_linkage.link_3.link_3_ground_reference_frame_ori &    
    j_marker_name = .four_bar_linkage.ground.ground_link_3_reference_frame_ori
! ------------------------------------------------------- !     
! ------------------------------------------------------- !