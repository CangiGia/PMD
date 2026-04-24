! +++++ CRANKSHAFT ::
! ------------------------------------------------------- !
! ------------------------------------------------------- !
part create &
	rigid_body name_and_position part_name=.crank_slider.crankshaft

part modify &
	rigid_body mass_properties part_name=.crank_slider.crankshaft material=.materials.steel

part attributes &
	part_name=.crank_slider.crankshaft color=YELLOW

marker create marker=.crank_slider.crankshaft.crankshaft_ground_reference_frame &
    location=0.0, 0.0, 0.0 &
    orientation=45.0, 0.0, 0.0

marker create marker=.crank_slider.crankshaft.crankshaft_rod_reference_frame &
    location=(LOC_RELATIVE_TO({(10.0cm),0.0,0.0}, .crank_slider.crankshaft.crankshaft_ground_reference_frame)) &
    orientation=(ORI_RELATIVE_TO({0.0,0.0,0.0}, .crank_slider.crankshaft.crankshaft_ground_reference_frame))

geometry create shape link &
    link_name=.crank_slider.crankshaft.crankshaft_geometry &
    width=(4.0cm) &
    depth=(2.0cm) &
    i_marker=.crank_slider.crankshaft.crankshaft_ground_reference_frame &
    j_marker=.crank_slider.crankshaft.crankshaft_rod_reference_frame
! ------------------------------------------------------- !
! ------------------------------------------------------- !

! +++++ ROD ::
! ------------------------------------------------------- !
! ------------------------------------------------------- !
part create &
	rigid_body name_and_position part_name=.crank_slider.rod

part modify &
	rigid_body mass_properties part_name=.crank_slider.rod material=.materials.steel

part attributes &
	part_name=.crank_slider.rod color=GREEN

marker create marker=.crank_slider.rod.rod_crankshaft_reference_frame &
    location=(LOC_RELATIVE_TO({0.0,0.0,0.0},.crank_slider.crankshaft.crankshaft_rod_reference_frame)) &
    orientation=(ORI_RELATIVE_TO({-(45+(90-79.8179)),0.0,0.0},.crank_slider.crankshaft.crankshaft_rod_reference_frame))

marker create marker=.crank_slider.rod.rod_slider_reference_frame &
    location=(LOC_RELATIVE_TO({(40.0cm),0.0,0.0}, .crank_slider.rod.rod_crankshaft_reference_frame)) &
    orientation=(ORI_RELATIVE_TO({0.0,0.0,0.0}, .crank_slider.rod.rod_crankshaft_reference_frame))

geometry create shape link &
    link_name=.crank_slider.rod.rod_geometry &
    width=(4.0cm) &
    depth=(2.0cm) &
    i_marker=.crank_slider.rod.rod_crankshaft_reference_frame &
    j_marker=.crank_slider.rod.rod_slider_reference_frame
! ------------------------------------------------------- !
! ------------------------------------------------------- !

! +++++ SLIDER ::
! ------------------------------------------------------- !
! ------------------------------------------------------- !
part create &
	rigid_body name_and_position part_name=.crank_slider.slider	

part modify &
	rigid_body mass_properties part_name=.crank_slider.slider material=.materials.steel

part attributes &
	part_name=.crank_slider.slider color=CYAN

marker create marker=.crank_slider.slider.slider_rod_reference_frame &
    location=(LOC_RELATIVE_TO({0.0,0.0,0.0},.crank_slider.rod.rod_slider_reference_frame)) &
    orientation=(ORI_RELATIVE_TO({0.0,0.0,0.0},.crank_slider.rod.rod_slider_reference_frame))

geometry create shape ellipsoid &
    ellipsoid_name=.crank_slider.slider.slider_geometry &
    x_scale_factor=(2*(2cm)) &
    y_scale_factor=(2*(2cm)) &
    z_scale_factor=(2*(2cm)) &
    center_marker=.crank_slider.slider.slider_rod_reference_frame
! ------------------------------------------------------- !
! ------------------------------------------------------- !