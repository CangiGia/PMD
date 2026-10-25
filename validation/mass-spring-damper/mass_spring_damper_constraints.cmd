! +++++ PRISMATIC JOINT :: GROUND - SLIDER
! ------------------------------------------------------- !
! ------------------------------------------------------- !
marker create marker=.mass_spring_damper.mass.mass_spring_damper_ground_ori &
    location=(LOC_RELATIVE_TO({0.0,0.0,0.0},.mass_spring_damper.mass.cm)) &
    orientation=(ORI_RELATIVE_TO({90,90,0.0},.mass_spring_damper.mass.cm))

marker create marker=.mass_spring_damper.ground.ground_mass_spring_damper &
    location=(LOC_RELATIVE_TO({0.0,0.0,0.0},.mass_spring_damper.mass.mass_spring_damper_ground_ori)) &
    orientation=(ORI_RELATIVE_TO({0.0,0.0,0.0},.mass_spring_damper.mass.mass_spring_damper_ground_ori))

constraint create joint translational &
    joint_name = .mass_spring_damper.mass_ground_joint &
    i_marker_name = .mass_spring_damper.mass.mass_spring_damper_ground_ori &    
    j_marker_name = .mass_spring_damper.ground.ground_mass_spring_damper

constraint modify joint translational joint_name=.mass_spring_damper.mass_ground_joint   &
      friction_enabled = "no" &
      translational_ic=200
! ------------------------------------------------------- !     
! ------------------------------------------------------- !