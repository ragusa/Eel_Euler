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

#include "ElementIntegralAbsVariablePostprocessor.h"

template<>
InputParameters validParams<ElementIntegralAbsVariablePostprocessor>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredParam<VariableName>("variable", "The name of the variable that this object operates on");
  return params;
}

ElementIntegralAbsVariablePostprocessor::ElementIntegralAbsVariablePostprocessor(const std::string & name, InputParameters parameters) :
    ElementIntegralPostprocessor(name, parameters),
    MooseVariableInterface(parameters, false),
    _var(_subproblem.getVariable(_tid, parameters.get<VariableName>("variable"))),
    _u(_var.sln()),
    _grad_u(_var.gradSln()),
    _u_dot(_var.uDot())
{
  addMooseVariableDependency(mooseVariable());
}

Real
ElementIntegralAbsVariablePostprocessor::computeQpIntegral()
{
    return std::fabs(_u[_qp]);
}
