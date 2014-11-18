
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
isJumpOn = false
Ce = 1.

###### Initial conditions ######
p_bc = 101325.
T_bc = 300.
gamma_bc = 0.
p0_bc = 101332.092927
T0_bc = 300.006

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
    type = StiffenedGasEquationOfState
  	gamma = 1.4
  	Pinf = 0.
  	q = 0.
  	Cv = 717.6
  	q_prime = 0. # reference entropy
  [../]

  [./JumpGradPress]
    type = JumpGradientInterface
    variable = pressure_aux
    jump_name = jump_grad_press_aux
  [../]

[]

###### Mesh #######
[Mesh]
#uniform_refine = 1
file = hump-2d.e
block_id = '1'
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
        value = 1.17666
	[../]
  [../]

  [./rhouA]
    family = LAGRANGE
    scaling = 1e-4
	[./InitialCondition]
        type = ConstantIC
        value = 4.0855247
	[../]
  [../]

  [./rhovA]
    family = LAGRANGE
    scaling = 1e-4
    [./InitialCondition]
    type = ConstantIC
    value = 0.
    [../]
   [../]

  [./rhoEA]
    family = LAGRANGE
    scaling = 1e-4
	[./InitialCondition]
        type = ConstantIC
        value = 253319.592751
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
    value = 0.05
    [../]
   [../]

   [./density_aux]
      family = LAGRANGE
	[./InitialCondition]
	type = ConstantIC
    value = 1.1766653
	[../]
   [../]

   [./total_energy_aux]
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
    value = 101325.
	[../]
   [../]

   [./temperature_aux]
    family = LAGRANGE
    [./InitialCondition]
    type = ConstantIC
    value = 300.
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

  [./TotEnerAK]
    type = TotalEnergyAux
    variable = total_energy_aux
    rhoEA = rhoEA
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

  [./TempAK]
    type = TemperatureAux
    variable = temperature_aux
    pressure = pressure_aux
    density = density_aux
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
    jump_grad_press = jump_grad_press_aux
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
  [../]

  [./mass_leak]
    type = NodalMassConservationPPS
    boundary = 'top bottom'
    rho_u = rhouA
    rho_v = rhovA
 [../]
[]

##############################################################################################
#                               BOUNDARY CONDITIONS                                          #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################
[NodalNormals]
    boundary = 'bottom top'
    corner_boundary = '100'
[]

[BCs]
  #active = ' '
  [./ContInflowDBC]
    type = EelStagnationPandTBC
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

  [./ContOutflowDBC]
    type = EelStaticPandTBC
    variable = rhoA
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    area = area_aux
    eos = eos
    equation_name = CONTINUITY
    boundary = 'right'
  [../]

  [./XMomInflowDBC]
    type = EelStagnationPandTBC
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

  [./XMomOutflowDBC]
    type = EelStaticPandTBC
    variable = rhouA
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    area = area_aux
    eos = eos
    equation_name = XMOMENTUM
    boundary = 'right'
  [../]

  [./XMomWallBC]
    type = MomentumFreeSlipBC
    variable = rhouA
    rho_u = rhouA
    rho_v = rhovA
    boundary = 'top bottom'
  [../]

  [./YMomInflowDBC]
    type = EelStagnationPandTBC
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

  [./YMomOutflowDBC]
    type = EelStaticPandTBC
    variable = rhovA
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    area = area_aux
    eos = eos
    equation_name = YMOMENTUM
    boundary = 'right'
  [../]

  [./YMomWallBC]
    type = MomentumFreeSlipBC
    variable = rhovA
    rho_u = rhouA
    rho_v = rhovA
    boundary = 'top bottom'
  [../]

  [./EnergyInflowDBC]
    type = EelStagnationPandTBC
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

  [./EnergyOutflowDBC]
    type = EelStaticPandTBC
    variable = rhoEA
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    rhoEA = rhoEA
    area = area_aux
    eos = eos
    equation_name = ENERGY
    boundary = 'right'
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
    solve_type = 'PJFNK'
#    petsc_options = '-snes_mf_operator -snes_ksp_ew'
#    petsc_options_iname = '-mat_fd_coloring_err  -mat_fd_type  -mat_mffd_type'
#    petsc_options_value = '1.e-12       ds             ds'
  [../]

  [./SMP_Newton]
    type = SMP
    full = true
    solve_type = 'PJFNK'
    petsc_options_iname = 'pc_type -pc_hypre_type'
    petsc_options_value = 'hypre boomerang'
  [../]
[]

##############################################################################################
#                                     EXECUTIONER                                            #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Executioner]
  type = Transient
  string scheme = 'implicit-euler' # ''bdf2'
  #num_steps = 20
  end_time = 2.
#  dt = 1.e-4
[./TimeStepper]
    type = FunctionDT
    time_t =  '0      1.e-2   2.e-2  1.'
    time_dt = '2.e-5  2.e-5   2.e-5  2.e-5'
  [../]
  dtmin = 1e-9
  #dtmax = 1e-5
  l_tol = 1e-8
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-5
  l_max_its = 50
  nl_max_its = 8
  [./Quadrature]
    type = GAUSS
    order = THIRD
  [../]
[]

##############################################################################################
#                                        OUTPUT                                              #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Output]
    output_initial = true
    file_base = Hump2DMach001EVs
    postprocessor_screen = false
    output_variables = 'pressure_aux density_aux norm_vel_aux kappa_aux kappa_max_aux mu_aux mu_max_aux mach_number_aux'
    interval = 200
    exodus = true
    #perf_log = true
[]
