! +++++ CREATE ROTATIONAL MOTION ::
! ------------------------------------------------------- !
! ------------------------------------------------------- !
constraint create motion motion_name=.four_bar_linkage.MOTION_1 &
    joint=.four_bar_linkage.ground_link_1_joint &
    type=rotational &
    time_derivative=displacement &
    function="25.0d * time"
! ------------------------------------------------------- !
! ------------------------------------------------------- !