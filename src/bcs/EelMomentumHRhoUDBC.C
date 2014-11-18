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

#include "EelMomentumHRhoUDBC.h"

template<>
InputParameters validParams<EelMomentumHRhoUDBC>()
{
  InputParameters params = validParams<NodalBC>();
  params.addRequiredParam<Real>("rhou", "Value of the momentum: rho*u");
  params.addRequiredCoupledVar("area", "area of the pipe");
  return params;
}

EelMomentumHRhoUDBC::EelMomentumHRhoUDBC(const std::string & name, InputParameters parameters) :
    NodalBC(name, parameters),
    _rhou(getParam<Real>("rhou")),
    _area(coupledValue("area"))
{}

Real
EelMomentumHRhoUDBC::computeQpResidual()
{
  return _u[_qp] - _area[_qp] * _rhou;
}
