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
This function computes the VariableTimeArea of the fluid.
**/
#include "VariableTimesAreaAux.h"

template<>
InputParameters validParams<VariableTimesAreaAux>()
{
  InputParameters params = validParams<AuxKernel>();
    params.addRequiredCoupledVar("var", "variable to multiply with area.");
    params.addRequiredCoupledVar("area", "area");
  return params;
}

VariableTimesAreaAux::VariableTimesAreaAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
    _var(coupledValue("var")),
    _area(coupledValue("area"))
{}

Real
VariableTimesAreaAux::computeValue()
{
  return _var[_qp] * _area[_qp];
}
