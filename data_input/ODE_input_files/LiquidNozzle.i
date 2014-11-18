#############################################################################
#                          USER OBJECTS                                     #
#############################################################################
# Define the user object class that store the EOS parameters.               #
#############################################################################

[UserObjects]
  [./eos]
    type = EquationOfState
  	gamma = 2.35
  	Pinf = 1.e9
  	q = 0.
  	Cv = 1816
  	q_prime = 0 # reference entropy
  [../]

[]

###### Mesh #######
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 100
  xmin = 0
  xmax = 1
  block_id = '0'
[]

###### Function #######
[Functions]
    [./ic]
        type = ParsedFunction
        value = 0.5+x*0.5
    [../]
[]

#############################################################################
#                             VARIABLES                                     #
#############################################################################
# Define the variables we want to solve for: l=liquid phase and g=gas phase.#
#############################################################################

[Variables]
  [./mach]
    family = LAGRANGE
    scaling = 1e+0
	[./InitialCondition]
        type = FunctionIC
        function = ic
	[../]
  [../]
[]

############################################################################################################
#                                            KERNELS                                                       #
############################################################################################################
# Define the kernels for time dependent, convection and viscosity terms. Same index as for variable block. #
############################################################################################################

[Kernels]
  [./MachK]
    type = EelFannoFlow
    variable = mach
    friction = 0.02
    Dh = 0.01
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
    type = DirichletBC
    variable = mach
    value = 0.5
    boundary = 'left'
  [../]

  [./ContOutflowDBC]
    type = DirichletBC
    variable = mach
    value = 1
    boundary = 'right'
  [../]
[]

##############################################################################################
#                                     EXECUTIONER                                            #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Executioner]
    type = Steady
    petsc_options = '-snes_mf_operator'
[]

##############################################################################################
#                                        OUTPUT                                              #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Output]
  output_initial = true
  #file_base =
  interval = 1
  exodus = true
  perf_log = true
[]
