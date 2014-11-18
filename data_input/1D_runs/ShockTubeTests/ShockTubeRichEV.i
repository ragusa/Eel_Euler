#
#####################################################
# Define some global parameters used in the blocks. #
#####################################################
#
[GlobalParams]
###### Boundary conditions: inflow and outflow #######
p0_bc = 1013250.0
T0_bc = 288.16
p_bc = 101325.0
T_bc = 288.16

###### Other parameters #######
order = FIRST
viscosity_name = ENTROPY
diffusion_name = ENTROPY
isJumpOn = true
Ce = 1.0

#Hw_fn = Hw_fn

###### Initial Conditions #######
pressure_init_left  = 1013250.
pressure_init_right = 101325.
vel_init_left = 0
vel_init_right = 0
temp_init_left = 288.16
temp_init_right = 288.16
membrane = 40
[]

##############################################################################################
#                                       FUNCTIONs                                            #
##############################################################################################
# Define functions that are used in the kernels and aux. kernels.                            #
##############################################################################################

[Functions]
#  [./Hw_fn]
#    type = ParsedFunction
#    value = 0.
#  [../]

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
    Pinf = 0
    q = 0.
    Cv = 717.645
    q_prime =  -23e2 # reference entropy
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

  [./JumpGradPressSmooth]
    type = SmoothFunction
    variable = jump_grad_press_aux
    var_name = jump_grad_press_smooth_aux
    execute_on = timestep_begin
  [../]

  [./JumpGradDensSmooth]
    type = SmoothFunction
    variable = jump_grad_dens_aux
    var_name = jump_grad_dens_smooth_aux
    execute_on = timestep_begin
  [../]

[]

###### Mesh #######
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 200
  xmin = 0
  xmax = 100
  block_id = '0'
#elem_type = EDGE3
[]

#############################################################################
#                             VARIABLES                                     #
#############################################################################
# Define the variables we want to solve for: l=liquid phase and g=gas phase.#
#############################################################################

[Variables]
  [./rhoA]
    family = LAGRANGE
    scaling = 1e-4
	[./InitialCondition]
        type = ConservativeVariables1DXIC
        eos = eos
        area = area
	[../]
  [../]

  [./rhouA]
    family = LAGRANGE
    scaling = 1e-4
    [./InitialCondition]
        type = ConservativeVariables1DXIC
        eos = eos
        area = area
    [../]
  [../]

  [./rhoEA]
    family = LAGRANGE
    scaling = 1e-6
	[./InitialCondition]
        type = ConservativeVariables1DXIC
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

  [./MomTime]
    type = EelTimeDerivative
    variable = rhouA
  [../]

  [./EnerTime]
    type = EelTimeDerivative
    variable = rhoEA
  [../]

  [./Mass]
    type = EelMass
    variable = rhoA
    rhouA_x = rhouA
  [../]

  [./Momentum]
    type = EelMomentum
    variable = rhouA
    rhoA = rhoA
    rhouA_x = rhouA
    rhoEA = rhoEA
    pressure = pressure_aux
    area = area_aux
    eos = eos
  [../]

  [./Energy]
    type = EelEnergy
    variable = rhoEA
    rhoA = rhoA
    rhouA_x = rhouA
    pressure = pressure_aux
    area = area_aux
    eos = eos
  [../]

  [./MassVisc]
    type = EelArtificialVisc
    variable = rhoA
    equation_name = CONTINUITY
    density = density_aux
    velocity_x = velocity_aux
    internal_energy = internal_energy_aux
    norm_velocity = norm_vel_aux
    area = area_aux
  [../]

   [./MomentumVisc]
    type = EelArtificialVisc
    variable = rhouA
    equation_name = XMOMENTUM
    density = density_aux
    velocity_x = velocity_aux
    internal_energy = internal_energy_aux
    norm_velocity = norm_vel_aux
    area = area_aux
  [../]

   [./EnergyVisc]
    type = EelArtificialVisc
    variable = rhoEA
    equation_name = ENERGY 
    density = density_aux
    velocity_x = velocity_aux
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

   [./velocity_aux]
      family = LAGRANGE
   [../]

   [./density_aux]
      family = LAGRANGE
   [../]

   [./internal_energy_aux]
      family = LAGRANGE
   [../]

   [./pressure_aux]
      family = LAGRANGE
   [../]

   [./mach_number_aux]
       family = LAGRANGE
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

  [./jump_grad_press_smooth_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./jump_grad_dens_smooth_aux]
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

  [./VelAK]
    type = VelocityAux
    variable = velocity_aux
    rhoA = rhoA
    rhouA = rhouA
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
    rhoEA = rhoEA
    area = area_aux
  [../]

  [./PressAK]
    type = PressureAux
    variable = pressure_aux
    rhoA = rhoA
    rhouA_x = rhouA
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
    area = area_aux
    eos = eos
  [../]

  [./NormVelAK]
    type = NormVectorAux
    variable = norm_vel_aux
    x_component = velocity_aux
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
    velocity_x = velocity_aux
    pressure = pressure_aux
    density = density_aux
    norm_velocity = norm_vel_aux
    jump_grad_press = jump_grad_press_smooth_aux
    jump_grad_dens = jump_grad_dens_smooth_aux
    eos = eos
    velocity_PPS_name = AverageVelocity
  [../]

[]

##############################################################################################
#                                     PPS                                                    #
##############################################################################################
# Define functions that are used in the kernels and aux. kernels.                            #
##############################################################################################
[Postprocessors]
#  [./MaxVelocity]
#    type = NodalMaxValue
#    variable = norm_vel_aux
#    execute_on = timestep
#  [../]

  [./AverageVelocity]
    type = ElementAverageValue
    variable = norm_vel_aux
#    execute_on = timestep
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
    value = 12.25
    boundary = 'left'
  [../]

  [./ContOutflowDBC]
    type = DirichletBC
    variable = rhoA
    value = 1.22
    boundary = 'right'
  [../]

  [./MomInflowDBC]
    type = DirichletBC
    variable = rhouA
    value = 0.
    boundary = 'left'
  [../]

  [./MomOutflowDBC]
    type = DirichletBC
    variable = rhouA
    value = 0.
    boundary = 'right'
  [../]

  [./EnergyInflowDBC]
    type = DirichletBC
    variable = rhoEA
    value = 2533125.
    boundary = 'left'
  [../]

  [./EnergyOutflowDBC]
    type = DirichletBC
    variable = rhoEA
    value = 253312.5
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
    petsc_options = '-snes_mf_operator -snes_ksp_ew'
    petsc_options_iname = '-mat_fd_coloring_err  -mat_fd_type  -mat_mffd_type'
    petsc_options_value = '1.e-12       ds             ds'
  [../]

  [./SMP_Newton]
    type = SMP
    full = true
    solve_type = 'PJFNK' # PJFNK, JFNK, NEWTON, FD
    line_search = 'none'
    #petsc_options = '-snes_mf'
    #petsc_options_iname = 'pc_type'
    #petsc_options_iname =  'pc_type -sub_pc_type'
    #petsc_options_value = 'lu'
    #petsc_options_value = 'ilu bjacobi'
  [../]
[]

##############################################################################################
#                                     EXECUTIONER                                            #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Executioner]
  type = Transient
  scheme = 'bdf2' # 'implicit-rk2'
  #rk_scheme = 'sdirk33'
  #num_steps = 1
  end_time = 0.08
  dt = 6.e-4
  [./TimeStepper]
    type = FunctionDT
    time_t =  '0      2.e-4  0.08'
    time_dt = '1.e-4  2.e-4  2.e-4'
  [../]
  dtmin = 1e-9
  l_tol = 1e-8
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-5
  l_max_its = 50
  nl_max_its = 10
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
#file_base = ToroTest1_sdirk22minus
  postprocessor_screen = false
  interval = 5
  exodus = true
  perf_log = true
[]
