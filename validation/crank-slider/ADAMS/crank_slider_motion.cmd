! +++++ CREATE ROTATIONAL MOTION ::
! ------------------------------------------------------- ! 
! ------------------------------------------------------- !
constraint create motion motion_name=.crank_slider.MOTION_1 &
    joint=.crank_slider.crankshaft_rod_joint &
    type=rotational &
    time_derivative=displacement &
    function="20.0d * time"
! ------------------------------------------------------- ! 
! ------------------------------------------------------- !