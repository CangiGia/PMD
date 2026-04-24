! +++++ REV JOINT :: GROUND - CRANKSHAFT
! ------------------------------------------------------- !
! ------------------------------------------------------- !
marker create marker=.crank_slider.ground.ground_crankshaft_reference_frame_ori &
    location=0.0, 0.0, 0.0 &
    orientation=0.0, 0.0, 0.0

marker create marker=.crank_slider.crankshaft.crankshaft_ground_reference_frame_ori &
    location=(LOC_RELATIVE_TO({0.0,0.0,0.0},.crank_slider.ground.ground_crankshaft_reference_frame_ori)) &
    orientation=(ORI_RELATIVE_TO({0.0,0.0,0.0},.crank_slider.ground.ground_crankshaft_reference_frame_ori))

constraint create joint revolute  &
	joint_name = .crank_slider.crankshaft_ground_joint  &
	friction_enabled = no  &
	i_marker_name = .crank_slider.ground.ground_crankshaft_reference_frame_ori &
	j_marker_name = .crank_slider.crankshaft.crankshaft_ground_reference_frame_ori
! ------------------------------------------------------- !
! ------------------------------------------------------- !

! +++++ REV JOINT :: CRANKSHAFT - ROD
! ------------------------------------------------------- ! 
! ------------------------------------------------------- !
marker create marker=.crank_slider.crankshaft.crankshaft_rod_reference_frame_ori &
    location=(LOC_RELATIVE_TO({0.0,0.0,0.0},.crank_slider.crankshaft.crankshaft_rod_reference_frame)) &
    orientation=(ORI_RELATIVE_TO({0.0,0.0,0.0},.crank_slider.crankshaft.crankshaft_ground_reference_frame_ori))

marker create marker=.crank_slider.rod.rod_crankshaft_reference_frame_ori &
    location=(LOC_RELATIVE_TO({0.0,0.0,0.0},.crank_slider.crankshaft.crankshaft_rod_reference_frame_ori)) &
    orientation=(ORI_RELATIVE_TO({0.0,0.0,0.0},.crank_slider.crankshaft.crankshaft_rod_reference_frame_ori))   

constraint create joint revolute  &
    joint_name = .crank_slider.crankshaft_rod_joint  &
    friction_enabled = no  &
    i_marker_name = .crank_slider.crankshaft.crankshaft_rod_reference_frame_ori &
    j_marker_name = .crank_slider.rod.rod_crankshaft_reference_frame_ori
! ------------------------------------------------------- !
! ------------------------------------------------------- !

! +++++ REV JOINT :: ROD - SLIDER
! ------------------------------------------------------- ! 
! ------------------------------------------------------- !
marker create marker=.crank_slider.rod.rod_slider_reference_frame_ori &
    location=(LOC_RELATIVE_TO({0.0,0.0,0.0},.crank_slider.slider.slider_rod_reference_frame)) &
    orientation=(ORI_RELATIVE_TO({0.0,0.0,0.0},.crank_slider.slider.slider_rod_reference_frame))

constraint create joint revolute  &
    joint_name = .crank_slider.rod_slider_joint  &
    friction_enabled = no  &
    i_marker_name = .crank_slider.slider.slider_rod_reference_frame &
    j_marker_name = .crank_slider.rod.rod_slider_reference_frame_ori
! ------------------------------------------------------- !
! ------------------------------------------------------- !

! +++++ PRISMATIC JOINT :: GROUND - SLIDER
! ------------------------------------------------------- !
! ------------------------------------------------------- !
marker create marker=.crank_slider.slider.slider_ground_reference_frame_ori &
    location=(LOC_RELATIVE_TO({0.0,0.0,0.0},.crank_slider.slider.slider_rod_reference_frame)) &
    orientation=(ORI_RELATIVE_TO({90,-90,0.0},.crank_slider.crankshaft.crankshaft_ground_reference_frame_ori))

marker create marker=.crank_slider.ground.ground_slider_reference_frame_ori &
    location=(LOC_RELATIVE_TO({0.0,0.0,0.0},.crank_slider.slider.slider_ground_reference_frame_ori)) &
    orientation=(ORI_RELATIVE_TO({0.0,0.0,0.0},.crank_slider.slider.slider_ground_reference_frame_ori))

constraint create joint translational &
    joint_name = .crank_slider.slider_ground_joint &
    i_marker_name = .crank_slider.slider.slider_ground_reference_frame_ori &    
    j_marker_name = .crank_slider.ground.ground_slider_reference_frame_ori
! ------------------------------------------------------- !     
! ------------------------------------------------------- !