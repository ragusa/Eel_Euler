[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 1
  nx = 10
 #elem_type = QUAD4
[]

[Variables]

  # ODE variables
  [./mach]
    family = SCALAR
    order = FIRST
    initial_condition = 1
  [../]

[]

[ScalarKernels]
  [./machSK]
    type = RayleighFannoFlow
    variable = mach
    friction = 0.02
    Dh = 0.01
  [../]
[]

[BCs]
  [./left]
    type = ScalarDirichletBC
    variable = mach
    boundary = 'left'
    scalar_var = 0.5
  [../]
[]

[Executioner]
  type = Transient
  start_time = 0
  dt = 0.01
  num_steps = 10
  petsc_options = '-snes_mf_operator'
[]

[Output]
  file_base = out
  output_initial = true
  exodus = true
[]
