! +++++ LINK-1 ::
! ------------------------------------------------------- !
! ------------------------------------------------------- !
! CREATE LINK
part create &
	rigid_body name_and_position part_name=.double_compound_pendulum.link_1

part modify &
	rigid_body mass_properties part_name=.double_compound_pendulum.link_1 material=.materials.steel

part attributes &
	part_name=.double_compound_pendulum.link_1 color=GREEN

marker create marker=.double_compound_pendulum.link_1.link_1_ground_reference_frame &
    location=0.0, 0.0, 0.0 &
    orientation=347.0053832081, 0.0, 0.0

marker create marker=.double_compound_pendulum.link_1.link_1_link_2_reference_frame &
    location=(LOC_RELATIVE_TO({(25.0cm),0.0,0.0}, .double_compound_pendulum.link_1.link_1_ground_reference_frame)) &
    orientation=347.0053832081, 0.0, 0.0

geometry create shape link &
    link_name=.double_compound_pendulum.link_1.link_1_geometry &
    width=(4.0cm) &
    depth=(2.0cm) &
    i_marker=.double_compound_pendulum.link_1.link_1_ground_reference_frame &
    j_marker=.double_compound_pendulum.link_1.link_1_link_2_reference_frame
! ------------------------------------------------------- !
! ------------------------------------------------------- !

! +++++ LINK-2 ::
! ------------------------------------------------------- !
! ------------------------------------------------------- !
part create &
	rigid_body name_and_position part_name=.double_compound_pendulum.link_2

part modify &
	rigid_body mass_properties part_name=.double_compound_pendulum.link_2 material=.materials.steel

part attributes &
	part_name=.double_compound_pendulum.link_2 color=YELLOW

marker create marker=.double_compound_pendulum.link_2.link_2_link_1_reference_frame &
    location=(LOC_RELATIVE_TO({0.0,0.0,0.0},.double_compound_pendulum.link_1.link_1_link_2_reference_frame)) &
	orientation=13.7697807428, 0.0, 0.0

marker create marker=.double_compound_pendulum.link_2.link_2_end_reference_frame &
    location=(LOC_RELATIVE_TO({(40.0cm),0.0,0.0},.double_compound_pendulum.link_2.link_2_link_1_reference_frame)) &
    orientation=(ORI_RELATIVE_TO({0.0,0.0,0.0},.double_compound_pendulum.link_2.link_2_link_1_reference_frame))

geometry create shape link &
    link_name=.double_compound_pendulum.link_2.link_2_geometry &
    width=(4.0cm) &
    depth=(2.0cm) &
    i_marker=.double_compound_pendulum.link_2.link_2_link_1_reference_frame &
    j_marker=.double_compound_pendulum.link_2.link_2_end_reference_frame
! ------------------------------------------------------- !
! ------------------------------------------------------- !