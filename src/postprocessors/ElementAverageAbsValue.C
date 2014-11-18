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

#include "ElementAverageAbsValue.h"

template<>
InputParameters validParams<ElementAverageAbsValue>()
{
  InputParameters params = validParams<ElementIntegralAbsVariablePostprocessor>();
  return params;
}

ElementAverageAbsValue::ElementAverageAbsValue(const std::string & name, InputParameters parameters) :
    ElementIntegralAbsVariablePostprocessor(name, parameters),
    _volume(0)
{}

void
ElementAverageAbsValue::initialize()
{
  ElementIntegralAbsVariablePostprocessor::initialize();
  _volume = 0;
}

void
ElementAverageAbsValue::execute()
{
  ElementIntegralAbsVariablePostprocessor::execute();

  _volume += _current_elem_volume;
}

Real
ElementAverageAbsValue::getValue()
{
  Real integral = ElementIntegralAbsVariablePostprocessor::getValue();

  gatherSum(_volume);

  return integral / _volume;
}

void
ElementAverageAbsValue::threadJoin(const UserObject & y)
{
  ElementIntegralAbsVariablePostprocessor::threadJoin(y);
  const ElementAverageAbsValue & pps = static_cast<const ElementAverageAbsValue &>(y);
  _volume += pps._volume;
}
