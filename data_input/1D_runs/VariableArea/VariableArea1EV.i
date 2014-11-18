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
Cmax = 1.
isLowMachShock = false

###### Initial Conditions #######
pressure_init_left = 400000.0
pressure_init_right = 200000.0
vel_init_left = 0.
vel_init_right = 329.130629
temp_init_left = 387.092778397
temp_init_right = 1221.60527572
membrane = 1.
length = 0.0
[]

##############################################################################################
#                                       FUNCTIONs                                            #
##############################################################################################
# Define functions that are used in the kernels and aux. kernels.                            #
##############################################################################################

[Functions]
  [./area]
    type = PiecewiseConstant
    x = '0.0 1.0 2.0'
    y = '1.0 0.6 0.6'
    axis='0 1'
  [../]

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
  	Cv = 717.6
  	q_prime = 0.# reference entropy
  [../]

  [./JumpGradArea]
    type = JumpGradientInterface
    variable = area_aux
    jump_name = jump_grad_area_aux
#    execute_on = timestep_begin
  [../]

  [./JumpGradPress]
    type = JumpGradientInterface
    variable = pressure_aux
    jump_name = jump_grad_press_aux
#    execute_on = timestep_begin
  [../]

  [./JumpGradDens]
    type = JumpGradientInterface
    variable = density_aux
    jump_name = jump_grad_dens_aux
#    execute_on = timestep_begin
  [../]

  [./JumpGradPressSmooth]
    type = SmoothFunction
    variable = jump_grad_press_aux
    var_name = jump_grad_press_smooth_aux
#    execute_on = timestep_begin
  [../]

  [./JumpGradDensSmooth]
    type = SmoothFunction
    variable = jump_grad_dens_aux
    var_name = jump_grad_dens_smooth_aux
#    execute_on = timestep_begin
  [../]
[]

###### Mesh #######
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 400
  xmin = 0.
  xmax = 2.
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
    scaling = 1e+0
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

  [./jump_grad_area_aux]
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
    type = FunctionAux
    variable = area_aux
    function = area
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
    jump_grad_area = jump_grad_area_aux
    area = area_aux
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
    value = 3.6
    boundary = 'left'
  [../]

  [./ContOutflowDBC]
    type = DirichletBC
    variable = rhoA
    value = 0.342222
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
    value = 112.635742118
    boundary = 'right'
  [../]

  [./EnergyInflowDBC]
    type = DirichletBC
    variable = rhoEA
    value = 1000000.
    boundary = 'left'
  [../]

  [./EnergyOutflowDBC]
    type = DirichletBC
    variable = rhoEA
    value = 318535.936326
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
    line_search = 'default'
  [../]
[]

##############################################################################################
#                                     EXECUTIONER                                            #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Executioner]
  type = Transient   # Here we use the Transient Executioner
  string scheme = 'bdf2'
  #rk_scheme = 'sdirk33'
  num_steps = 1000
  end_time = 7.e-4
  dt = 1.e-5
  [./TimeStepper]
    type = FunctionDT
    time_t =  '0      1.e-6  0.2'
    time_dt = '1.e-6  1.e-6  1.e-6'
  [../]
  dtmin = 1e-9
  l_tol = 1e-8
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-6
  l_max_its = 50
  nl_max_its = 10
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
    interval = 1
    console = true
    exodus = true
[]
