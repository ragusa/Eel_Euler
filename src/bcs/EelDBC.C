/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "EelDBC.h"

template<>
InputParameters validParams<EelDBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  return params;
}

EelDBC::EelDBC(const std::string & name, InputParameters parameters) :
    IntegratedBC(name, parameters)
{}

Real
EelDBC::computeQpResidual()
{
    return -_grad_u[_qp]*_normals[_qp]*_test[_i][_qp];
}

Real
EelDBC::computeQpJacobian()
{
  // TODO
  return 0;
}

Real
EelDBC::computeQpOffDiagJacobian(unsigned jvar)
{
  // TODO
  return 0;
}
