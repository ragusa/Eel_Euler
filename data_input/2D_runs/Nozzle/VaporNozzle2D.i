#
#####################################################
# Define some global parameters used in the blocks. #
#####################################################
#

[GlobalParams]
###### Boundary conditions: inflow and outflow #######
p0_bc = 1e6
T0_bc = 453.
p_bc = 0.5e6
T_bc = 453.

###### Other parameters #######
order = FIRST
viscosity_name = FIRST_ORDER
diffusion_name = ENTROPY 

###### Initial Conditions #######
pressure_init_left = 1e6
pressure_init_right = 0.5e6
vel_init_left = 0
vel_init_right = 0
temp_init_left = 453
temp_init_right = 453
membrane = 0.5
[]

#############################################################################
#                          USER OBJECTS                                     #
#############################################################################
# Define the user object class that store the EOS parameters.               #
#############################################################################

[UserObjects]
  [./eos]
    type = EquationOfState
  	gamma = 1.34
  	Pinf = 0.
  	q = 1968e3
  	Cv = 1265
  	q_prime = 0. # reference entropy
  [../]

[]

###### Mesh #######
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 50
  ny = 1
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
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
    scaling = 1e-1
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
        type = ConstantIC
        value = 0.
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
    vel_x = velocity_x_aux
    vel_y = velocity_y_aux
    pressure = pressure_aux
    area = area_aux
    component = 0
  [../]

  [./YMomentum]
    type = EelMomentum
    variable = rhovA
    vel_x = velocity_x_aux
    vel_y = velocity_y_aux
    pressure = pressure_aux
    area = area_aux
    component = 1
  [../]

  [./Energy]
    type = EelEnergy
    variable = rhoEA
    rhoA = rhoA
    rhouA_x = rhouA
    rhouA_y = rhovA
    pressure = pressure_aux
    area = area_aux
  [../]

  [./MassVisc]
    type = EelArtificialVisc
    variable = rhoA
    equation_name = CONTINUITY
    density = density_aux
    velocity_x = vel_area_x_aux
    velocity_y = vel_area_y_aux
    internal_energy = internal_energy_aux
    norm_velocity = norm_vel_aux
    area = area_aux
  [../]

   [./XMomentumVisc]
    type = EelArtificialVisc
    variable = rhouA
    equation_name = XMOMENTUM
    density = density_aux
    velocity_x = vel_area_x_aux
    velocity_y = vel_area_y_aux
    internal_energy = internal_energy_aux
    norm_velocity = norm_vel_aux
    area = area_aux
  [../]

  [./YMomentumVisc]
    type = EelArtificialVisc
    variable = rhovA
    equation_name = YMOMENTUM
    density = density_aux
    velocity_x = vel_area_x_aux
    velocity_y = vel_area_y_aux
    internal_energy = internal_energy_aux
    norm_velocity = norm_vel_aux
    area = area_aux
  [../]

   [./EnergyVisc]
    type = EelArtificialVisc
    variable = rhoEA
    equation_name = ENERGY 
    density = density_aux
    velocity_x = vel_area_x_aux
    velocity_y = vel_area_y_aux
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
   [../]

   [./velocity_y_aux]
    family = LAGRANGE
   [../]

   [./vel_area_x_aux]
    family = LAGRANGE
   [../]

   [./vel_area_y_aux]
    family = LAGRANGE
   [../]

   [./mach_number_aux]
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

  [./VelAreaXAK]
    type = VariableTimesAreaAux
    variable = vel_area_x_aux
    var = velocity_x_aux
    area = area_aux
  [../]

  [./VelAreaYAK]
    type = VariableTimesAreaAux
    variable = vel_area_y_aux
    var = velocity_y_aux
    area = area_aux
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

  [./TemperatureAK]
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
    block = '0'
    velocity_x = velocity_x_aux
    velocity_y = velocity_y_aux
    pressure = pressure_aux
    density = density_aux
    norm_velocity = norm_vel_aux
    eos = eos
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
    area = area_aux
    eos = eos
    #type = DirichletBC
    #value = 7.698810226
    boundary = 'left'
  [../]

  [./ContOutflowDBC]
    type = EelStaticPandTBC
    variable = rhoA
    equation_name = CONTINUITY
    vel_x = velocity_x_aux
    vel_y = velocity_y_aux
    temperature = temperature_aux
    area = area_aux
    eos = eos
    #    type = DirichletBC
    #    value = 3.849405113
    boundary = 'right'
  [../]

  [./ContTopBC]
    type = EelWallBC
    variable = rhoA
    pressure = pressure_aux
#    velocity_x = velocity_x_aux
#    velocity_y = velocity_y_aux
    area = area_aux
    equation_name = CONTINUITY
    boundary = 'top'
  [../]

  [./ContBottomflowBC]
    type = EelWallBC
    variable = rhoA
    pressure = pressure_aux
#    velocity_x = velocity_x_aux
#    velocity_y = velocity_y_aux
    area = area_aux
    equation_name = CONTINUITY
    boundary = 'bottom'
  [../]

  [./XMomInflowDBC]
    type = EelStagnationPandTBC
    variable = rhouA
    equation_name = XMOMENTUM
    rhoA = rhoA
    rhouA_x = rhouA
    area = area_aux
    eos = eos
    #    type = DirichletBC
    #    value = 0
    boundary = 'left'
  [../]

  [./XMomOutflowDBC]
    type = EelStaticPandTBC
    variable = rhouA
    equation_name = XMOMENTUM
    vel_x = velocity_x_aux
    vel_y = velocity_y_aux
    temperature = temperature_aux
    area = area_aux
    eos = eos
    #   type = DirichletBC
    #   value = 0
    boundary = 'right'
  [../]

  [./XMomTopBC]
    type = EelWallBC
    variable = rhouA
    pressure = pressure_aux
#    velocity_x = velocity_x_aux
#    velocity_y = velocity_y_aux
    area = area_aux
    equation_name = XMOMENTUM
    boundary = 'top'
  [../]

  [./XMomBottomBC]
    type = EelWallBC
    variable = rhouA
    pressure = pressure_aux
#    velocity_x = velocity_x_aux
#    velocity_y = velocity_y_aux
    area = area_aux
    equation_name = XMOMENTUM
    boundary = 'bottom'
  [../]

  [./YMomInflowDBC]
    type = EelStagnationPandTBC
    variable = rhovA
    equation_name = YMOMENTUM
    rhoA = rhoA
    rhouA_x = rhouA
    area = area_aux
    eos = eos
    #    type = DirichletBC
    #    value = 0
    boundary = 'left'
  [../]

  [./YMomOutflowDBC]
    type = EelStaticPandTBC
    variable = rhovA
    equation_name = YMOMENTUM
    vel_x = velocity_x_aux
    vel_y = velocity_y_aux
    temperature = temperature_aux
    area = area_aux
    eos = eos
    boundary = 'right'
  [../]

  [./YMomTopBC]
    type = EelWallBC
    variable = rhovA
    pressure = pressure_aux
#    velocity_x = velocity_x_aux
#    velocity_y = velocity_y_aux
    area = area_aux
    equation_name = YMOMENTUM
    boundary = 'top'
  [../]

  [./YMomBottomBC]
    type = EelWallBC
    variable = rhovA
    pressure = pressure_aux
#    velocity_x = velocity_x_aux
#    velocity_y = velocity_y_aux
    area = area_aux
    equation_name = YMOMENTUM
    boundary = 'bottom'
  [../]

  [./EnergyInflowDBC]
    type = EelStagnationPandTBC
    variable = rhoEA
    equation_name = ENERGY
    rhoA = rhoA
    rhouA_x = rhouA
    area = area_aux
    eos = eos
    boundary = 'left'
  [../]

  [./EnergyOutflowDBC]
    type = EelStaticPandTBC
    variable = rhoEA
    equation_name = ENERGY
    vel_x = velocity_x_aux
    vel_y = velocity_y_aux
    temperature = temperature_aux
    area = area_aux
    eos = eos
    boundary = 'right'
  [../]

  [./EnergyTopBC]
    type = EelWallBC
    variable = rhoEA
    pressure = pressure_aux
#    velocity_x = velocity_x_aux
#    velocity_y = velocity_y_aux
    area = area_aux
    equation_name = ENERGY
    boundary = 'top'
  [../]

  [./EnergyBottomBC]
    type = EelWallBC
    variable = rhoEA
    pressure = pressure_aux
#    velocity_x = velocity_x_aux
#    velocity_y = velocity_y_aux
    area = area_aux
    equation_name = ENERGY 
    boundary = 'bottom'
  [../]

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
     Ao = 1.0
     Bo = 0.0
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
    petsc_options = '-snes_mf_operator -snes_ksp_ew'
    petsc_options_iname = '-mat_fd_coloring_err'
    petsc_options_value = '1.e-12'
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
  type = Transient   # Here we use the Transient Executioner
  #scheme = 'explicit-euler'
  string scheme = 'bdf2'
  #petsc_options = '-snes'
  #petsc_options_iname = '-pc_type'
  #petsc_options_value = 'lu'
  #num_steps = 100
  end_time = 2.6e-2
  dt = 1e-4
  dtmin = 1e-9
  #dtmax = 1e-5
  l_tol = 1e-8
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-5
  l_max_its = 50
  nl_max_its = 60
  #predictor_scale = 0.0
  #e_tol = 0.01
  #e_max = 0.05
[]

##############################################################################################
#                                        OUTPUT                                              #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Output]
  output_initial = true
  interval = 1
  exodus = true
  perf_log = true
[]
