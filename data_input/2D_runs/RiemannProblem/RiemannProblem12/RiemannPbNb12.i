# Website with Initial conditions for 2D Riemann problems:
# http://www.astrophysik.uni-kiel.de/index.php?option=com_content&view=article&id=262&Itemid=168&lang=en
#####################################################
# Define some global parameters used in the blocks. #
#####################################################
#
[GlobalParams]
###### Other parameters #######
order = FIRST
viscosity_name = ENTROPY
diffusion_name = PARABOLIC # ENTROPY
isJumpOn = true
Ce = 1.

###### Initial conditions ######
press_upper_left_corner = 1.
press_upper_right_corner = 0.4
press_bottom_left_corner = 1.
press_bottom_right_corner = 1.

dens_upper_left_corner = 1.
dens_upper_right_corner = 0.53125
dens_bottom_left_corner = 0.8
dens_bottom_right_corner = 1.

x_vel_upper_left_corner = 0.7276068751
x_vel_upper_right_corner = 0.
x_vel_bottom_left_corner = 0.
x_vel_bottom_right_corner = 0.

y_vel_upper_left_corner = 0.
y_vel_upper_right_corner = 0.
y_vel_bottom_left_corner = 0.
y_vel_bottom_right_corner = 0.7276068751

x_node = 0.5
y_node = 0.5
[]

##############################################################################################
#                                       FUNCTIONs                                            #
##############################################################################################
# Define functions that are used in the kernels and aux. kernels.                            #
##############################################################################################

[Functions]
  [./area]
    type = ParsedFunction
    value = 1.
  [../]
[]

#############################################################################
#                          USER OBJECTS                                     #
#############################################################################
# Define the user object class that store the EOS parameters.               #
#############################################################################

[UserObjects]
  [./eos]
    type = StiffenedGasEquationOfState
  	gamma = 1.4
  	Pinf = 0.
  	q = 0.
  	Cv = 2.5
  	q_prime = 0. # reference entropy
  [../]

  [./JumpGradPress]
    type = JumpGradientInterface
    variable = pressure_aux
    jump_name = jump_grad_press_aux
    execute_on = timestep_begin
  [../]

  [./JumpGradDens]
    type = JumpGradientInterface
    variable = density_aux
    jump_name = jump_grad_dens_aux
    execute_on = timestep_begin
  [../]

  [./SmoothJumpGradPress]
    type = SmoothFunction
    variable = jump_grad_press_aux
    var_name = smooth_jump_grad_press_aux
    execute_on = timestep_begin
  [../]

  [./SmoothJumpGradDens]
    type = SmoothFunction
    variable = jump_grad_dens_aux
    var_name = smooth_jump_grad_dens_aux
    execute_on = timestep_begin
  [../]

[]

###### Mesh #######
[Mesh]
    type = GeneratedMesh
    uniform_refine = 0
    dim = 2
    xmin = 0.
    xmax = 1
    ymin = 0.
    ymax = 1.
    nx = 100
    ny = 100
[]

#############################################################################
#                             VARIABLES                                     #
#############################################################################
# Define the variables we want to solve for: l=liquid phase and g=gas phase.#
#############################################################################

[Variables]
  [./rhoA]
    family = LAGRANGE
    scaling = 1e+1
	[./InitialCondition]
        type = FourSquaresIC2D
        eos = eos
        area = area
	[../]
  [../]

  [./rhouA]
    family = LAGRANGE
    scaling = 1e+1
    [./InitialCondition]
        type = FourSquaresIC2D
        eos = eos
        area = area
    [../]
  [../]

  [./rhovA]
    family = LAGRANGE
    scaling = 1e+1
    [./InitialCondition]
        type = FourSquaresIC2D
        eos = eos
        area = area
    [../]
   [../]

  [./rhoEA]
    family = LAGRANGE
    scaling = 1e+1
    [./InitialCondition]
        type = FourSquaresIC2D
        eos = eos
        area = area
    [../]
  [../]
[]

############################################################################################################
#                                            KERNELS                                                       #
############################################################################################################
# Define the kernels for time dependent, convection and viscosity terms. Same index as for variable block. #
############################################################################################################

[Kernels]

  [./ContTime]
    type = EelTimeDerivative
    variable = rhoA
  [../]

  [./XMomTime]
    type = EelTimeDerivative
    variable = rhouA
  [../]

  [./YMomTime]
    type = EelTimeDerivative
    variable = rhovA
  [../]

  [./EnerTime]
    type = EelTimeDerivative
    variable = rhoEA
  [../]

  [./Mass]
    type = EelMass
    variable = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
  [../]

  [./XMomentum]
    type = EelMomentum
    variable = rhouA
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    pressure = pressure_aux
    area = area_aux
    component = 0
    eos = eos
  [../]

  [./YMomentum]
    type = EelMomentum
    variable = rhovA
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    pressure = pressure_aux
    area = area_aux
    component = 1
    eos = eos
  [../]

  [./Energy]
    type = EelEnergy
    variable = rhoEA
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    pressure = pressure_aux
    area = area_aux
    eos = eos
  [../]

  [./MassVisc]
    type = EelArtificialVisc
    variable = rhoA
    equation_name = CONTINUITY
    density = density_aux
    velocity_x = velocity_x_aux
    velocity_y = velocity_y_aux
    internal_energy = internal_energy_aux
    norm_velocity = norm_vel_aux
    area = area_aux
  [../]

   [./XMomentumVisc]
    type = EelArtificialVisc
    variable = rhouA
    equation_name = XMOMENTUM
    density = density_aux
    velocity_x = velocity_x_aux
    velocity_y = velocity_y_aux
    internal_energy = internal_energy_aux
    norm_velocity = norm_vel_aux
    area = area_aux
  [../]

  [./YMomentumVisc]
    type = EelArtificialVisc
    variable = rhovA
    equation_name = YMOMENTUM
    density = density_aux
    velocity_x = velocity_x_aux
    velocity_y = velocity_y_aux
    internal_energy = internal_energy_aux
    norm_velocity = norm_vel_aux
    area = area_aux
  [../]

   [./EnergyVisc]
    type = EelArtificialVisc
    variable = rhoEA
    equation_name = ENERGY 
    density = density_aux
    velocity_x = velocity_x_aux
    velocity_y = velocity_y_aux
    internal_energy = internal_energy_aux
    norm_velocity = norm_vel_aux
    area = area_aux
  [../]
[]

##############################################################################################
#                                       AUXILARY VARIABLES                                   #
##############################################################################################
# Define the auxilary variables                                                              #
##############################################################################################

[AuxVariables]

   [./area_aux]
      family = LAGRANGE
   [../]

   [./velocity_x_aux]
      family = LAGRANGE
	[./InitialCondition]
	type = ConstantIC
    value = 0.
	[../]
   [../]

   [./velocity_y_aux]
    family = LAGRANGE
    [./InitialCondition]
    type = ConstantIC
    value = 0.
    [../]
   [../]

   [./mach_number_aux]
    family = LAGRANGE
    [./InitialCondition]
    type = ConstantIC
    value = 0.
    [../]
   [../]

   [./density_aux]
      family = LAGRANGE
	[./InitialCondition]
	type = ConstantIC
    value = 0.
	[../]
   [../]

   [./internal_energy_aux]
      family = LAGRANGE
	[./InitialCondition]
	type = ConstantIC
    value = 0.
	[../]
   [../]

   [./pressure_aux]
      family = LAGRANGE
	[./InitialCondition]
	type = ConstantIC
    value = 0.5e6
	[../]
   [../]

   [./norm_vel_aux]
    family = LAGRANGE
   [../]

   [./mu_max_aux]
    family = MONOMIAL
    order = CONSTANT
   [../]

   [./kappa_max_aux]
    family = MONOMIAL
    order = CONSTANT
   [../]

   [./mu_aux]
    family = MONOMIAL
    order = CONSTANT
   [../]

   [./kappa_aux]
    family = MONOMIAL
    order = CONSTANT
   [../]

  [./jump_grad_press_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./jump_grad_dens_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./smooth_jump_grad_press_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./smooth_jump_grad_dens_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]
##############################################################################################
#                                       AUXILARY KERNELS                                     #
##############################################################################################
# Define the auxilary kernels for liquid and gas phases. Same index as for variable block.   #
##############################################################################################

[AuxKernels]

  [./AreaAK]
    type = AreaAux
    variable = area_aux
    area = area
  [../]

  [./VelXAK]
    type = VelocityAux
    variable = velocity_x_aux
    rhoA = rhoA
    rhouA = rhouA
  [../]

  [./VelYAK]
    type = VelocityAux
    variable = velocity_y_aux
    rhoA = rhoA
    rhouA = rhovA
  [../]

  [./DensAK]
    type = DensityAux
    variable = density_aux
    rhoA = rhoA
    area = area_aux
  [../]

  [./IntEnerAK]
    type = InternalEnergyAux
    variable = internal_energy_aux
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    area = area_aux
  [../]

  [./PressAK]
    type = PressureAux
    variable = pressure_aux
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    area = area_aux
    eos = eos
  [../]

  [./MachNumAK]
    type = MachNumberAux
    variable = mach_number_aux
    pressure = pressure_aux
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    area = area_aux
    eos = eos
  [../]

  [./NormVelAK]
    type = NormVectorAux
    variable = norm_vel_aux
    x_component = velocity_x_aux
    y_component = velocity_y_aux
   [../]

   [./MuMaxAK]
    type = MaterialRealAux
    variable = mu_max_aux
    property = mu_max
   [../]

   [./KappaMaxAK]
    type = MaterialRealAux
    variable = kappa_max_aux
    property = kappa_max
   [../]

   [./MuAK]
    type = MaterialRealAux
    variable = mu_aux
    property = mu
   [../]

   [./KappaAK]
    type = MaterialRealAux
    variable = kappa_aux
    property = kappa
   [../]
[]
##############################################################################################
#                                       MATERIALS                                            #
##############################################################################################
# Define functions that are used in the kernels and aux. kernels.                            #
##############################################################################################

[Materials]
#active = ''
  [./EntViscMat]
    type = ComputeViscCoeff
    block = '0'
    velocity_x = velocity_x_aux
    velocity_y = velocity_y_aux
    pressure = pressure_aux
    density = density_aux
    norm_velocity = norm_vel_aux
    jump_grad_press = smooth_jump_grad_press_aux
    jump_grad_dens = smooth_jump_grad_dens_aux
    velocity_PPS_name = MaxVelocity
    eos = eos
  [../]
[]
##############################################################################################
#                                     PPS                                                    #
##############################################################################################
# Define functions that are used in the kernels and aux. kernels.                            #
##############################################################################################
[Postprocessors]
  [./MaxVelocity]
    type = NodalMaxValue
    variable = norm_vel_aux
    execute_on = timestep_begin
  [../]
[]
##############################################################################################
#                               BOUNDARY CONDITIONS                                          #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################
[BCs]
  #active = ' '
  [./ContLeftBC]
    type = EelFluxBC
    variable = rhoA
    equation_name = CONTINUITY
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    area = area_aux
    eos = eos
    boundary = 'left'
  [../]

  [./ContRightBC]
    type = EelFluxBC
    variable = rhoA
    equation_name = CONTINUITY
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    area = area_aux
    eos = eos
    boundary = 'right'
  [../]

  [./ContTopBC]
    type = EelFluxBC
    variable = rhoA
    equation_name = CONTINUITY
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    area = area_aux
    eos = eos
    boundary = 'top'
  [../]

  [./ContBottomBC]
    type = EelFluxBC
    variable = rhoA
    equation_name = CONTINUITY
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    area = area_aux
    eos = eos
    boundary = 'bottom'
  [../]

  [./XMomLeftBC]
    type = EelFluxBC
    variable = rhouA
    equation_name = XMOMENTUM
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    area = area_aux
    eos = eos
    boundary = 'left'
  [../]

  [./XMomRightBC]
    type = EelFluxBC
    variable = rhouA
    equation_name = XMOMENTUM
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    area = area_aux
    eos = eos
    boundary = 'right'
  [../]

  [./XMomTopBC]
    type = EelFluxBC
    variable = rhouA
    equation_name = XMOMENTUM
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    area = area_aux
    eos = eos
    boundary = 'top'
  [../]

  [./XMomBottomBC]
    type = EelFluxBC
    variable = rhouA
    equation_name = XMOMENTUM
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    area = area_aux
    eos = eos
    boundary = 'bottom'
  [../]

  [./YMomLeftBC]
    type = EelFluxBC
    variable = rhovA
    equation_name = YMOMENTUM
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    area = area_aux
    eos = eos
    boundary = 'left'
  [../]

  [./YMomRightBC]
    type = EelFluxBC
    variable = rhovA
    equation_name = YMOMENTUM
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    area = area_aux
    eos = eos
    boundary = 'right'
  [../]

  [./YMomTopBC]
    type = EelFluxBC
    variable = rhovA
    equation_name = YMOMENTUM
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    area = area_aux
    eos = eos
    boundary = 'top'
  [../]

  [./YMomBottomBC]
    type = EelFluxBC
    variable = rhovA
    equation_name = YMOMENTUM
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    area = area_aux
    eos = eos
    boundary = 'bottom'
  [../]

  [./EnergyLeftBC]
    type = EelFluxBC
    variable = rhoEA
    equation_name = ENERGY
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    area = area_aux
    eos = eos
    boundary = 'left'
  [../]

  [./EnergyRightBC]
    type = EelFluxBC
    variable = rhoEA
    equation_name = ENERGY
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    area = area_aux
    eos = eos
    boundary = 'right'
  [../]

  [./EnergyTopBC]
    type = EelFluxBC
    variable = rhoEA
    equation_name = ENERGY
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    area = area_aux
    eos = eos
    boundary = 'top'
  [../]

  [./EnergyBottomBC]
    type = EelFluxBC
    variable = rhoEA
    equation_name = ENERGY
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    area = area_aux
    eos = eos
    boundary = 'bottom'
  [../]
[]

##############################################################################################
#                                  PRECONDITIONER                                            #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Preconditioning]
#active = 'FDP_Newton'
    active = 'SMP_Newton'
  [./FDP_Newton]
    type = FDP
    full = true
    petsc_options = '-snes_mf_operator -snes_ksp_ew'
    petsc_options_iname = '-mat_fd_coloring_err  -mat_fd_type  -mat_mffd_type'
    petsc_options_value = '1.e-12       ds             ds'
  [../]

  [./SMP_Newton]
    type = SMP
    full = true
    solve_type = 'PJFNK'
    line_search = 'none'
  [../]
[]

##############################################################################################
#                                     EXECUTIONER                                            #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Executioner]
  type = Transient
  scheme = 'bdf2'
#num_steps = 20
  end_time = 0.2
    #dt = 1.e-4
[./TimeStepper]
    type = FunctionDT
    time_t =  '0      2.e-6  0.2'
    time_dt = '1.e-6  1.e-3  1.e-3'
  [../]
  dtmin = 1e-9
  #dtmax = 1e-5
  l_tol = 1e-8
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-5
  l_max_its = 50
  nl_max_its = 10
  [./Quadrature]
    type = TRAP
  [../]

#[Adaptivity]
#  initial_marker = errorfrac
#  initial_steps = 3
#  marker = errorfrac
#  max_h_level = 5
#  [./Indicators]
#    [./error]
##        type = GradientJumpIndicator
#        type = FluxJumpIndicator
#        #type = LaplacianJumpIndicator
#        variable = rhoA
#        property = kappa
#        block = '0'
#    [../]
#  [../]
#  [./Markers]
#    [./errorfrac]
#        type = ErrorFractionMarker
#        #type = ErrorToleranceMarker
#        refine = 0.3
#        coarsen = 0.1
#        indicator = error
#    [../]
#  [../]
# [../]
[]

##############################################################################################
#                                        OUTPUT                                              #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Output]
  output_initial = true
#file_base = CompressionCorner2D
  postprocessor_screen = false
  interval = 1
  exodus = true
  perf_log = true
[]
