! +++++ CREATE TRANSLATIONAL SPRING DAMPER ::
! ------------------------------------------------------- ! 
! ------------------------------------------------------- !	
marker create marker=.mass_spring_damper.ground.mass_spring_damper_force &
    location=0.0, 0.0, 0.0 &
    orientation=0.0, 0.0, 0.0

marker create marker=.mass_spring_damper.mass.mass_spring_damper_force &
    location=(LOC_RELATIVE_TO({0.0,800,0.0},.mass_spring_damper.ground.mass_spring_damper_force)) &
    orientation=(ORI_RELATIVE_TO({0.0,0.0,0.0},.mass_spring_damper.ground.mass_spring_damper_force))
! ------------------------------------------------------- ! 
! ------------------------------------------------------- !

undo begin
default coordinate_system default_coordinate_system=.mass_spring_damper
undo begin suppress=yes
   assembly create instance instance=.mass_spring_damper.spring &
      definition=.mdi.forces.spring
   variable modify variable=.mass_spring_damper.spring.i_marker obj=(.mass_spring_damper.mass.mass_spring_damper_force)
   variable modify variable=.mass_spring_damper.spring.j_marker obj=(.mass_spring_damper.ground.mass_spring_damper_force)
   variable modify variable=.mass_spring_damper.spring.stiffness_mode str="linear"
   variable modify variable=.mass_spring_damper.spring.stiffness_coefficient real=(4.0E-03(newton/mm))
   variable modify variable=.mass_spring_damper.spring.stiffness_spline obj=(none)
   variable modify variable=.mass_spring_damper.spring.damping_mode str="linear"
   variable modify variable=.mass_spring_damper.spring.damping_coefficient real=(1.0E-03(newton-sec/mm))
   variable modify variable=.mass_spring_damper.spring.damping_spline obj=(none)
   variable modify variable=.mass_spring_damper.spring.free_length_mode str="design_length"
   variable modify variable=.mass_spring_damper.spring.i_dynamic_visibility str="on"
   variable modify variable=.mass_spring_damper.spring.j_dynamic_visibility str="off"
   entity attributes entity_name=.mass_spring_damper.spring color=red
undo end
group modify group=SELECT_LIST object=.mass_spring_damper.spring
undo end

var set var=.mass_spring_damper.spring.stiffness_mode          str="linear"
var set var=.mass_spring_damper.spring.stiffness_coefficient   real=(1.0E3(newton/mm))
var set var=.mass_spring_damper.spring.damping_mode            str="linear"
var set var=.mass_spring_damper.spring.damping_coefficient     real=(10(newton-sec/mm))
var set var=.mass_spring_damper.spring.free_length_mode        str="Constant_Value"
var set var=.mass_spring_damper.spring.free_length             real=(EVAL(DM(.mass_spring_damper.mass.mass_spring_damper_force, .mass_spring_damper.ground.mass_spring_damper_force)))
var set var=.mass_spring_damper.spring.preload                 real=0
var set var=.mass_spring_damper.spring.spring_visibility       str="depends"
var set var=.mass_spring_damper.spring.damper_visibility       str="depends"
var set var=.mass_spring_damper.spring.i_dynamic_visibility    str="On"
var set var=.mass_spring_damper.spring.j_dynamic_visibility    str="Off"
! -------------------------------------------------------------- !
! -------------------------------------------------------------- !