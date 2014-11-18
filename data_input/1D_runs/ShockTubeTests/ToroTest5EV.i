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
Cjump = 5.
isShock = true

###### Initial Conditions #######
pressure_init_left = 1000.
pressure_init_right = 0.01
vel_init_left = -19.59745
vel_init_right = -19.59745
temp_init_left = 1000.
temp_init_right = 0.01
membrane = 0.8
length = 0.
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
  	Pinf = 0
  	q = 0.
  	Cv = 2.5
  	q_prime = 0.# reference entropy
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
  nx = 500
  xmin = 0
  xmax = 1
  block_id = '0'
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
        type = ConservativeVariables1DXIC
        eos = eos
        area = area
	[../]
  [../]

  [./rhouA]
    family = LAGRANGE
    scaling = 1e-3
    [./InitialCondition]
        type = ConservativeVariables1DXIC
        eos = eos
        area = area
    [../]
  [../]

  [./rhoEA]
    family = LAGRANGE
    scaling = 1e-3
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
    jump_grad_press = jump_grad_press_aux
    jump_grad_dens = jump_grad_dens_aux
    eos = eos
    rhov2_PPS_name = AverageRhovel2
    rhoc2_PPS_name = AverageRhoc2
  [../]

[]

##############################################################################################
#                                     PPS                                                    #
##############################################################################################
# Define functions that are used in the kernels and aux. kernels.                            #
##############################################################################################
[Postprocessors]
[./AverageRhovel2]
    type = ElementAverageMultipleValues
    variable = norm_vel_aux
    output_type = RHOVEL2
    rhoA = rhoA
    rhouA_x = rhouA
    rhoEA = rhoEA
    eos = eos
    area = area_aux
[../]

[./AverageRhoc2]
    type = ElementAverageMultipleValues
    variable = norm_vel_aux
    output_type = RHOC2
    rhoA = rhoA
    rhouA_x = rhouA
    rhoEA = rhoEA
    eos = eos
    area = area_aux
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
    value = 1.
    boundary = 'left'
  [../]

  [./ContOutflowDBC]
    type = DirichletBC
    variable = rhoA
    value = 1.
    boundary = 'right'
  [../]

  [./MomInflowDBC]
    type = DirichletBC
    variable = rhouA
    value = -19.59745
    boundary = 'left'
  [../]

  [./MomOutflowDBC]
    type = DirichletBC
    variable = rhouA
    value = -19.59745
    boundary = 'right'
  [../]

  [./EnergyInflowDBC]
    type = DirichletBC
    variable = rhoEA
    value = 2692.030023
    boundary = 'left'
  [../]

  [./EnergyOutflowDBC]
    type = DirichletBC
    variable = rhoEA
    value = 192.0550233
    boundary = 'right'
  [../]
[]

##############################################################################################
#                                  PRECONDITIONER                                            #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Preconditioning]
 active = 'FDP_Newton'
#    active = 'SMP_Newton'
  [./FDP_Newton]
    type = FDP
    full = true
    solve_type = 'PJFNK'
#    petsc_options = '-snes_mf_operator -snes_ksp_ew'
    petsc_options_iname = '-mat_fd_coloring_err  -mat_fd_type  -mat_mffd_type'
    petsc_options_value = '1.e-12       ds             ds'
  [../]

  [./SMP_Newton]
    type = SMP
    full = true
    solve_type= 'PJFNK'
    line_search = 'none'
#    petsc_options = '-snes'
#    petsc_options_iname = 'pc_type'
#    petsc_options_value = 'ilu'
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
  #num_steps = 400
  end_time = 0.012
  #dt = 5e-6
  [./TimeStepper]
    type = FunctionDT
    time_t =  '0     2.e-6  0.012'
    time_dt = '1.e-6 2.e-5  2.e-5'
  [../]
  dtmin = 1e-9
  #dtmax = 1e-5
  l_tol = 1e-8
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  l_max_its = 50
  nl_max_its = 30
  [./Quadrature]
    type = TRAP
    order = SECOND
  [../]
[]
##############################################################################################
#                                        OUTPUT                                              #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Outputs]
    output_initial = true
    postprocessor_screen = false
    interval = 1
    console = true
    exodus = true
    perf_log = true
[]
