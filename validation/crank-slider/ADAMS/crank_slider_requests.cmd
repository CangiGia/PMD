! +++++ CONSTRAINT VIOLATION REQUESTS ::
! ------------------------------------------------------- !
! ------------------------------------------------------- !
output_control create request &
    request_name = .crank_slider.req_constraints_violation &
    f1 = "TM(.crank_slider.ground.ground_crankshaft_reference_frame_ori,.crank_slider.crankshaft.crankshaft_ground_reference_frame_ori)" &
    f2 = "TM(.crank_slider.crankshaft.crankshaft_rod_reference_frame_ori,.crank_slider.rod.rod_crankshaft_reference_frame_ori)" &
    f3 = "TM(.crank_slider.slider.slider_rod_reference_frame,.crank_slider.rod.rod_slider_reference_frame_ori)" &
    f4 = "TM(.crank_slider.slider.slider_ground_reference_frame_ori,.crank_slider.ground.ground_slider_reference_frame_ori)" &
    f5 = "0.0" &
    f6 = "0.0" &
    f7 = "0.0" &
    f8 = "0.0" 
! ------------------------------------------------------- !
! ------------------------------------------------------- !



! +++++ COG COORDINATES REQUEST ::
! ------------------------------------------------------- !
! ------------------------------------------------------- !
output_control create request &
    request_name = .crank_slider.req_cm_positions &
    f1 = "DX(.crank_slider.crankshaft.cm, .crank_slider.ground.ground_crankshaft_reference_frame_ori)" &
    f2 = "DY(.crank_slider.crankshaft.cm, .crank_slider.ground.ground_crankshaft_reference_frame_ori)" &
    f3 = "DX(.crank_slider.rod.cm, .crank_slider.ground.ground_crankshaft_reference_frame_ori)" &
    f4 = "DY(.crank_slider.rod.cm, .crank_slider.ground.ground_crankshaft_reference_frame_ori)" &
    f5 = "DX(.crank_slider.slider.cm, .crank_slider.ground.ground_crankshaft_reference_frame_ori)" &
    f6 = "DY(.crank_slider.slider.cm, .crank_slider.ground.ground_crankshaft_reference_frame_ori)" &
    f7 = "0.0" &
    f8 = "0.0" 
! ------------------------------------------------------- !
! ------------------------------------------------------- !
