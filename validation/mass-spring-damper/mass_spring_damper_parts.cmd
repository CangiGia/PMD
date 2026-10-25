! +++++ MASS ::
! ------------------------------------------------------- !
! ------------------------------------------------------- !
part create rigid_body name_and_position part_name=.mass_spring_damper.mass
part modify rigid_body mass_properties part_name=.mass_spring_damper.mass material=.materials.steel
part attributes part_name=.mass_spring_damper.mass color=RED name_vis=off

! marker create marker=.mass_spring_damper.mass.mass_spring_damper &
    ! location=-250.0, 500.0, 0.0 &
    ! orientation=0.0, 0.0, 0.0
	
marker create marker=.mass_spring_damper.mass.mass_spring_damper &
    location=-250, 800, -250 &
    orientation=0.0, 0.0, 0.0
	
geometry create shape block &
    block_name=.mass_spring_damper.mass.mass &
    diag_corner_coords=( 500mm ), ( 500mm ), ( 500mm ) &
    corner_marker=.mass_spring_damper.mass.mass_spring_damper

part modify rigid mass_properties  &
  part_name = .mass_spring_damper.mass  &
  mass = 10000  &
  center_of_mass_marker = .mass_spring_damper.mass.cm  &
  ixx = 0.001  &
  iyy = 0.001  &
  izz = 0.001  &
! ------------------------------------------------------- !
! ------------------------------------------------------- !