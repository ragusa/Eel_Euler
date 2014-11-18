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

#include "MassMatrix.h"

template<>
InputParameters validParams<MassMatrix>()
{
  InputParameters params = validParams<TimeDerivative>();
  return params;
}

MassMatrix::MassMatrix(const std::string & name,
                                             InputParameters parameters) :
    TimeDerivative(name,parameters)
{}

Real
MassMatrix::computeQpResidual()
{
    return _u[_qp] * _test[_i][_qp];
}

Real
MassMatrix::computeQpJacobian()
{
    return _test[_j][_qp] * _test[_i][_qp];
}
