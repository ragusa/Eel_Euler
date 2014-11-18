#
#####################################################
# Define some global parameters used in the blocks. #
#####################################################
#
[GlobalParams]
###### Other parameters #######
order = FIRST
viscosity_name = ENTROPY
diffusion_name = ENTROPY
isJumpOn = true
Ce = 1.

###### Initial conditions ######
p_bc = 5
T_bc = 1.4
gamma_bc = 0.

Hw_fn = Hw_fn
[]

##############################################################################################
#                                       FUNCTIONs                                            #
##############################################################################################
# Define functions that are used in the kernels and aux. kernels.                            #
##############################################################################################

[Functions]
  [./Hw_fn]
    type = ParsedFunction
    value = 0.
  [../]

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
    type = EquationOfState
  	gamma = 1.4
  	Pinf = 0.
  	q = 0.
  	Cv =  2.5
  	q_prime = 0. # reference entropy
  [../]

  [./JumpGradPress]
    type = JumpGradientInterface
    variable = pressure_aux
    jump_name = jump_grad_press_aux
  [../]

  [./SmoothJumpGradPress]
    type = SmoothFunction
    variable = jump_grad_press_aux
    var_name = smooth_jump_grad_press_aux
  [../]

  [./JumpGradDens]
    type = JumpGradientInterface
    variable = density_aux
    jump_name = jump_grad_dens_aux
  [../]

  [./SmoothJumpGradDens]
    type = SmoothFunction
    variable = jump_grad_dens_aux
    var_name = smooth_jump_grad_dens_aux
  [../]
[]

###### Mesh #######
[Mesh]
  uniform_refine = 4
  file = compression_corner.e
  block_id = '1'
  boundary_id = '1 2 3'
  boundary_name = 'wall outflow inflow'
[]

#############################################################################
#                             VARIABLES                                     #
#############################################################################
# Define the variables we want to solve for: l=liquid phase and g=gas phase.#
#############################################################################

[Variables]
  [./rhoA]
    family = LAGRANGE
    scaling = 1e+0
	[./InitialCondition]
        type = ConstantIC
        value = 7.
	[../]
  [../]

  [./rhouA]
    family = LAGRANGE
    scaling = 1e+0
	[./InitialCondition]
        type = ConstantIC
        value = 17.5
	[../]
  [../]

  [./rhovA]
    family = LAGRANGE
    scaling = 1e+0
    [./InitialCondition]
    type = ConstantIC
    value = 0.
    [../]
   [../]

  [./rhoEA]
    family = LAGRANGE
    scaling = 1e+0
	[./InitialCondition]
        type = ConstantIC
        value = 34.375
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
    block = '1'
    velocity_x = velocity_x_aux
    velocity_y = velocity_y_aux
    pressure = pressure_aux
    density = density_aux
    norm_velocity = norm_vel_aux
    jump_grad_press = smooth_jump_grad_press_aux
    jump_grad_press = smooth_jump_grad_dens_aux
    pressure_PPS_name = AveragePressure
    velocity_PPS_name = AverageVelocity
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
  [../]

  [./AverageVelocity]
    type = ElementAverageValue
    variable = norm_vel_aux
  [../]
[]

##############################################################################################
#                               BOUNDARY CONDITIONS                                          #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################
[BCs]
  #active = ' '
  [./ContInflowDBC]
    type = DirichletBC
    variable = rhoA
    value = 7.
    boundary = 'inflow'
  [../]

  [./ContOutflowDBC]
    type = EelStaticPandTBC
    variable = rhoA
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    density = density_aux
    pressure = pressure_aux
    vel_x = velocity_x_aux
    vel_y = velocity_y_aux
    area = area_aux
    eos = eos
    equation_name = CONTINUITY
    boundary = 'outflow'
  [../]

  [./ContWallBC]
    type = EelWallBC
    variable = rhoA
    pressure = pressure_aux
    area = area_aux
    eos = eos
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    equation_name = CONTINUITY
    boundary = 'wall'
  [../]

  [./XMomInflowDBC]
    type = DirichletBC
    variable = rhouA
    value = 17.5
    boundary = 'inflow'
  [../]

  [./XMomOutflowDBC]
    type = EelStaticPandTBC
    variable = rhouA
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    density = density_aux
    pressure = pressure_aux
    vel_x = velocity_x_aux
    vel_y = velocity_y_aux
    area = area_aux
    eos = eos
    equation_name = XMOMENTUM
    boundary = 'outflow'
  [../]

  [./XMomWallBC]
    type = EelWallBC
    variable = rhouA
    pressure = pressure_aux
    area = area_aux
    eos = eos
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    equation_name = XMOMENTUM
    boundary = 'wall'
  [../]

  [./YMomInflowDBC]
    type = DirichletBC
    variable = rhovA
    value = 0.
    boundary = 'inflow'
  [../]

  [./YMomOutflowDBC]
    type = EelStaticPandTBC
    variable = rhovA
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    density = density_aux
    pressure = pressure_aux
    vel_x = velocity_x_aux
    vel_y = velocity_y_aux
    area = area_aux
    eos = eos
    equation_name = YMOMENTUM
    boundary = 'outflow'
  [../]

  [./YMomWallBC]
    type = EelWallBC
    variable = rhovA
    pressure = pressure_aux
    area = area_aux
    eos = eos
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    equation_name = YMOMENTUM
    boundary = 'wall'
  [../]

  [./EnergyInflowDBC]
    type = DirichletBC
    variable = rhoEA
    value = 34.375
    boundary = 'inflow'
  [../]

  [./EnergyOutflowDBC]
    type = EelStaticPandTBC
    variable = rhoEA
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    density = density_aux
    pressure = pressure_aux
    vel_x = velocity_x_aux
    vel_y = velocity_y_aux
    area = area_aux
    eos = eos
    equation_name = ENERGY
    boundary = 'outflow'
  [../]

  [./EnergyWallBC]
    type = EelWallBC
    variable = rhoEA
    pressure = pressure_aux
    area = area_aux
    eos = eos
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    equation_name = ENERGY
    boundary = 'wall'
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
    #petsc_options = '-snes_mf_operator -ksp_converged_reason -ksp_monitor -snes_ksp_ew'
    #petsc_options_iname = '-pc_type'
    #petsc_options_value = 'lu'
  [../]

  [./SMP_Newton]
    type = SMP
    full = true
    solve_type = 'PJFNK'
    line_search = 'none'
#petsc_options = '-snes'
#petsc_options_iname = 'pc_type'
#petsc_options_value = 'lu'
  [../]
[]

##############################################################################################
#                                     EXECUTIONER                                            #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Executioner]
  type = Transient
  string scheme = 'bdf2'
  #num_steps = 100
  end_time = 1.5
  dt = 1.e-4
[./TimeStepper]
    type = FunctionDT
    time_t =  '0      2.e-6  1.5'
    time_dt = '1.e-6  1.e-3  1.e-3'
  [../]
  dtmin = 1e-9
  #dtmax = 1e-5
  l_tol = 1e-8
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-5
  l_max_its = 50
  nl_max_its = 8
  [./Quadrature]
    type = TRAP
  [../]

#[Adaptivity]
#    marker = errorfrac
#    max_h_level = 5
#    [./Indicators]
#        [./error]
#        type = GradientJumpIndicator
#        variable = rhoA
#       [../]
#    [../]
#    [./Markers]
#        [./errorfrac]
#        block = '1'
#        type = ErrorFractionMarker
#       refine = 0.5
#        coarsen = 0.1
#        indicator = error
#        [../]
#    [../]
#[../]

[]

##############################################################################################
#                                        OUTPUT                                              #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Output]
  output_initial = true
  file_base = CompressionCorner2D
  postprocessor_screen = false
  interval = 4
  exodus = true
  perf_log = true
[]
