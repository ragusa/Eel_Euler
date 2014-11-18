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
/**
This function computes the density of the fluid.
**/
#include "DensityAux.h"

template<>
InputParameters validParams<DensityAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("rhoA", "rhoA");
  params.addRequiredCoupledVar("area", "area");
  return params;
}

DensityAux::DensityAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
    _rhoA(coupledValue("rhoA")),
    _area(coupledValue("area"))
{}

Real
DensityAux::computeValue()
{
  return _rhoA[_qp] / _area[_qp];
}
