! +++++ LINK-1 ::
! ------------------------------------------------------- !
! ------------------------------------------------------- !
part create &
	rigid_body name_and_position part_name=.four_bar_linkage.link_1

part modify &
	rigid_body mass_properties part_name=.four_bar_linkage.link_1 material=.materials.steel

part attributes &
	part_name=.four_bar_linkage.link_1 color=GREEN

marker create marker=.four_bar_linkage.link_1.link_1_ground_reference_frame &
    location=0.0, 0.0, 0.0 &
    orientation=40.0, 0.0, 0.0

marker create marker=.four_bar_linkage.link_1.link_1_link_2_reference_frame &
    location=(LOC_RELATIVE_TO({(40.0mm),0.0,0.0},.four_bar_linkage.link_1.link_1_ground_reference_frame)) &
    orientation=(ORI_RELATIVE_TO({0.0,0.0,0.0},.four_bar_linkage.link_1.link_1_ground_reference_frame))

geometry create shape link &
    link_name=.four_bar_linkage.link_1.link_1_geometry &
    width=(8.0mm) &
    depth=(4.0mm) &
    i_marker=.four_bar_linkage.link_1.link_1_ground_reference_frame &
    j_marker=.four_bar_linkage.link_1.link_1_link_2_reference_frame
! ------------------------------------------------------- !
! ------------------------------------------------------- !

! +++++ LINK-2 ::
! ------------------------------------------------------- !
! ------------------------------------------------------- !
part create &
	rigid_body name_and_position part_name=.four_bar_linkage.link_2

part modify &
	rigid_body mass_properties part_name=.four_bar_linkage.link_2 material=.materials.steel

part attributes &
	part_name=.four_bar_linkage.link_2 color=YELLOW

marker create marker=.four_bar_linkage.link_2.link_2_link_1_reference_frame &
    location=(LOC_RELATIVE_TO({0.0,0.0,0.0},.four_bar_linkage.link_1.link_1_link_2_reference_frame)) &
    orientation=20.298, 0.0, 0.0

marker create marker=.four_bar_linkage.link_2.link_2_link_3_reference_frame &
    location=(LOC_RELATIVE_TO({(120.0mm),0.0,0.0},.four_bar_linkage.link_2.link_2_link_1_reference_frame)) &
    orientation=(ORI_RELATIVE_TO({0.0,0.0,0.0},.four_bar_linkage.link_2.link_2_link_1_reference_frame))

geometry create shape link &
    link_name=.four_bar_linkage.link_2.link_2_geometry &
    width=(8.0mm) &
    depth=(4.0mm) &
    i_marker=.four_bar_linkage.link_2.link_2_link_1_reference_frame &
    j_marker=.four_bar_linkage.link_2.link_2_link_3_reference_frame
! ------------------------------------------------------- !
! ------------------------------------------------------- !

! +++++ LINK-3 ::
! ------------------------------------------------------- !
! ------------------------------------------------------- !
part create &
	rigid_body name_and_position part_name=.four_bar_linkage.link_3

part modify &
	rigid_body mass_properties part_name=.four_bar_linkage.link_3 material=.materials.steel

part attributes &
	part_name=.four_bar_linkage.link_3 color=RED

marker create marker=.four_bar_linkage.link_3.link_3_link_2_reference_frame &
    location=(LOC_RELATIVE_TO({0.0,0.0,0.0},.four_bar_linkage.link_2.link_2_link_3_reference_frame)) &
    orientation=-122.605, 0.0, 0.0

marker create marker=.four_bar_linkage.link_3.link_3_ground_reference_frame &
    location=(LOC_RELATIVE_TO({(80.0mm),0.0,0.0},.four_bar_linkage.link_3.link_3_link_2_reference_frame)) &
    orientation=(ORI_RELATIVE_TO({0.0,0.0,0.0},.four_bar_linkage.link_3.link_3_link_2_reference_frame))

geometry create shape link &
    link_name=.four_bar_linkage.link_3.link_3_geometry &
    width=(8.0mm) &
    depth=(4.0mm) &
    i_marker=.four_bar_linkage.link_3.link_3_link_2_reference_frame &
    j_marker=.four_bar_linkage.link_3.link_3_ground_reference_frame
! ------------------------------------------------------- !
! ------------------------------------------------------- !