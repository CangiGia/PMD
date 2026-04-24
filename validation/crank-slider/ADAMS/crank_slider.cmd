! +++++ DEFAULT SCRIPT LAYOUT FOR MODEL CREATION :: 
! ------------------------------------------------------- !
! ------------------------------------------------------- !
! "" The following rows can be considered as standard. 
! Pay attention to change the name of the model during copy 
! and paste operations. ""
! ------------------------------------------------------- !
! ------------------------------------------------------- !



! +++++ DEFAULT :: 
! ------------------------------------------------------- !
! ------------------------------------------------------- !
! +++++ MODEL DECLARATION :: 
model create &
	model_name="crank_slider" &

! +++++ MODEL DEFAULT UNITS :: 
defaults units  &
	length = mm  &
	angle = deg  &
	force = newton &
	mass = kg    &
	time = sec

defaults units  &
	coordinate_system_type = cartesian  &
	orientation_type = body313

model attributes model = .crank_slider size_of_icons = 20

! +++++ GRAVITY :: 
force create body gravitational gravity = Gravity &
	x_comp = 0 &
	y_comp = -9806.65 &
	z_comp = 0

entity attributes &
	entity_name     = .crank_slider.Gravity &
	type_filter     = Gravity_Field &
	visibility      = no_opinion &
	name_visibility = no_opinion &
	entity_scope    = all_color &
	scale           = 2 &
	transparency    = 0

! +++++ GRID SETTING VIEW :: 
interface grid modify &
	dots_vis    = yes &
	dot_color   = No_Color &
	dot_size    = 1 &
	lines_vis   = no &
	line_weight = 1 &
	line_style  = Solid &
	line_color  = No_Color &
	axes_vis    = yes &
	axis_weight = 1 &
	axis_color  = No_Color &
	triad_vis   = no &
	extent      = ((25mm)),((25mm)) &
	spacing     = ((25mm)),((25mm)) &
! ------------------------------------------------------- !
! ------------------------------------------------------- !



! +++++ RIGID MODEL CREATION (macros calling) ::
! +++++ PARTS CREATION :: 
! ------------------------------------------------------- !
! ------------------------------------------------------- !
file command read file = "crank_slider_parts.cmd"
! ------------------------------------------------------- !
! ------------------------------------------------------- !



! +++++ CONSTRAINTS CREATION :: 
! ------------------------------------------------------- !
! ------------------------------------------------------- !
file command read file = "crank_slider_constraints.cmd"
! ------------------------------------------------------- !
! ------------------------------------------------------- !



! +++++ MOTION CREATION ::
! ------------------------------------------------------- !
! ------------------------------------------------------- !
file command read file = "crank_slider_motion.cmd"
! ------------------------------------------------------- !
! ------------------------------------------------------- !



! +++++ REQUESTS CREATION ::
! ------------------------------------------------------- !
! ------------------------------------------------------- !
file command read file = "crank_slider_requests.cmd"
! ------------------------------------------------------- !
! ------------------------------------------------------- !