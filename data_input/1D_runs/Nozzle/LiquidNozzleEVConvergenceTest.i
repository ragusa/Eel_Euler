#
#####################################################
# Define some global parameters used in the blocks. #
#####################################################
#

[GlobalParams]
###### Boundary conditions: inflow and outflow #######
p0_bc = 1.e6
T0_bc = 453.
p_bc = 0.5e6
T_bc = 453.

###### Other parameters #######
order = FIRST
viscosity_name = ENTROPY
diffusion_name = ENTROPY
isJumpOn = true
Ce = 1.
Cjump = 0.

###### Initial Conditions #######
pressure_init_left = 1.e6
pressure_init_right = 0.5e6
vel_init_left = 0
vel_init_right = 0
temp_init_left = 453
temp_init_right = 453
membrane = 0.5
length = 1.
[]

#############################################################################
#                          USER OBJECTS                                     #
#############################################################################
# Define the user object class that store the EOS parameters.               #
#############################################################################

[UserObjects]
  [./eos]
    type = StiffenedGasEquationOfState
  	gamma = 2.35
  	Pinf = 1.e9
  	q = -1167e3
  	Cv = 1816
  	q_prime = 0 # reference entropy
  [../]

  [./JumpGradPress]
    type = JumpGradientInterface
    variable = pressure_aux
    jump_name = jump_grad_press_aux
  [../]
[]

###### Mesh #######
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 16
  xmin = 0
  xmax = 1
  block_id = '0'
[]

##############################################################################################
#                                       FUNCTIONs                                            #
##############################################################################################
# Define functions that are used in the kernels and aux. kernels.                            #
##############################################################################################

[Functions]

  #active = 'area'
  [./area]
    type = AreaFunction
    #value = Ao * ( 1 + 0.5*cos((x-left)/l*pi) ) + Bo
    left = 0.0
    length = 1.
    Ao = 1.0
    Bo = 0.0
  [../]
  
  [./exact_sol_press]
    type = ExactSolAreaVariable
    variable_name = PRESSURE
    length = 1.
    eos = eos
  [../]

  [./exact_sol_dens]
    type = ExactSolAreaVariable
    variable_name = DENSITY
    length = 1.
    eos = eos
  [../]

  [./exact_sol_vel]
    type = ExactSolAreaVariable
    variable_name = VELOCITY
    length = 1.
    eos = eos
  [../]

[]

#############################################################################
#                             VARIABLES                                     #
#############################################################################
# Define the variables we want to solve for: l=liquid phase and g=gas phase.#
#############################################################################

[Variables]
  [./rhoA]
    family = LAGRANGE
    scaling = 1e-1
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
        type = ConstantIC
        value = 0.
	[../]
  [../]

  [./rhoEA]
    family = LAGRANGE
    scaling = 1e-4
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

   [./press_exact_aux]
      family = LAGRANGE
   [../]

   [./vel_exact_aux]
    family = LAGRANGE
   [../]

   [./dens_exact_aux]
    family = LAGRANGE
   [../]

   [./area_aux]
      family = LAGRANGE
   [../]

   [./velocity_aux]
      family = LAGRANGE
   [../]

   [./density_aux]
      family = LAGRANGE
   [../]

   [./total_energy_aux]
      family = LAGRANGE
   [../]

   [./internal_energy_aux]
      family = LAGRANGE
   [../]

   [./pressure_aux]
      family = LAGRANGE
   [../]

   [./temperature_aux]
    family = LAGRANGE
   [../]

   [./mach_number_aux]
      family = LAGRANGE
   [../]

   [./norm_vel_aux]
    family = LAGRANGE
   [../]

   [./vel_area_aux]
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

  [./ExactSolPressAK]
    type = FunctionAux
    variable = press_exact_aux
    function = exact_sol_press
  [../]

  [./ExactSolDensityAK]
    type = FunctionAux
    variable = dens_exact_aux
    function = exact_sol_dens
  [../]

  [./ExactSolVelAK]
    type = FunctionAux
    variable = vel_exact_aux
    function = exact_sol_vel
  [../]

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
    area = area_aux
    eos = eos
  [../]

  [./NormVelAK]
    type = NormVectorAux
    variable = norm_vel_aux
    x_component = velocity_aux
  [../]

  [./VariableTimesAreaAK]
    type = VariableTimesAreaAux
    variable = vel_area_aux
    var = velocity_aux
    area = area_aux
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

  [./L2ErrorPressure]
    type = ElementL2Error
    variable = pressure_aux
    function = exact_sol_press
  [../]

  [./L2ErrorDensity]
    type = ElementL2Error
    variable = density_aux
    function = exact_sol_dens
  [../]

  [./L2ErrorVel]
    type = ElementL2Error
    variable = velocity_aux
    function = exact_sol_vel
  [../]

  [./L1ErrorPressure]
    type = ElementL1Error
    variable = pressure_aux
    function = exact_sol_press
  [../]

  [./L1ErrorDensity]
    type = ElementL1Error
    variable = density_aux
    function = exact_sol_dens
  [../]

  [./L1ErrorVel]
    type = ElementL1Error
    variable = velocity_aux
    function = exact_sol_vel
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
    type = EelStagnationPandTBC
    variable = rhoA
    equation_name = CONTINUITY
    rhoA = rhoA
    rhouA_x = rhouA
    rhoEA = rhoEA
    area = area_aux
    eos = eos
    boundary = 'left'
  [../]

  [./ContOutflowDBC]
    type = EelStaticPandTBC
    variable = rhoA
    equation_name = CONTINUITY
    rhoA = rhoA
    rhouA_x = rhouA
    rhoEA = rhoEA
    vel_x = velocity_aux    
    pressure = pressure_aux
    density = density_aux
    area = area_aux
    eos = eos
    boundary = 'right'
  [../]

  [./MomInflowDBC]
    type = EelStagnationPandTBC
    variable = rhouA
    equation_name = XMOMENTUM
    rhoA = rhoA
    rhouA_x = rhouA
    rhoEA = rhoEA
    area = area_aux
    eos = eos
    boundary = 'left'
  [../]

  [./MomOutflowDBC]
    type = EelStaticPandTBC
    variable = rhouA
    equation_name = XMOMENTUM
    rhoA = rhoA
    rhouA_x = rhouA
    rhoEA = rhoEA
    vel_x = velocity_aux
    pressure = pressure_aux
    density = density_aux
    area = area_aux
    eos = eos
    boundary = 'right'
  [../]

  [./EnergyInflowDBC]
    type = EelStagnationPandTBC
    variable = rhoEA
    equation_name = ENERGY
    rhoA = rhoA
    rhouA_x = rhouA
    rhoEA = rhoEA
    area = area_aux
    eos = eos
    boundary = 'left'
  [../]

  [./EnergyOutflowDBC]
    type = EelStaticPandTBC
    variable = rhoEA
    equation_name = ENERGY
    rhoA = rhoA
    rhouA_x = rhouA
    rhoEA = rhoEA
    vel_x = velocity_aux
    pressure = pressure_aux
    density = density_aux
    area = area_aux
    eos = eos
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
  [./FDP_Newton]
    type = FDP
    full = true
    solve_type = NEWTON
    petsc_options = '-snes_mf_operator -snes_ksp_ew'
    petsc_options_iname = '-mat_fd_coloring_err  -mat_fd_type  -mat_mffd_type'
    petsc_options_value = '1.e-12       ds             ds'
    #petsc_options = '-snes_mf_operator -ksp_converged_reason -ksp_monitor -snes_ksp_ew'
    #petsc_options_iname = '-pc_type'
    #petsc_options_value = 'lu'
  [../]

  [./SMP]
  type=SMP
    full=true
    petsc_options = '-snes_mf_operator'
    petsc_options_iname = '-pc_type'
    petsc_options_value = 'lu'
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
  end_time = 1.
  dt = 1.e-3
#  [./TimeStepper]
#    type = FunctionDT
#    time_t =  '0     2.6e-2'
#    time_dt = '5e-5  5e-5'
#  [../]
  dtmin = 1e-9
  #dtmax = 1e-5
  l_tol = 1e-8
  nl_rel_tol = 1e-7
  nl_abs_tol = 1e-6
  l_max_its = 50
  nl_max_its = 30
  [./Quadrature]
    type = TRAP
    order = TENTH
  [../]
[]

##############################################################################################
#                                        OUTPUT                                              #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Outputs]
    #file_base = LiquidNozzleEV_nel_8
    output_initial = true
    interval = 1
    console = true
    exodus = true
    postprocessor_screen = false
    perf_log = true
[]
