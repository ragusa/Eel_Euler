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

#include "ElementAverageMultipleValues.h"

template<>
InputParameters validParams<ElementAverageMultipleValues>()
{
  InputParameters params = validParams<ElementIntegralMultipleVariablesPostprocessor>();
  return params;
}

ElementAverageMultipleValues::ElementAverageMultipleValues(const std::string & name, InputParameters parameters) :
    ElementIntegralMultipleVariablesPostprocessor(name, parameters),
    _volume(0)
{}

void
ElementAverageMultipleValues::initialize()
{
  ElementIntegralMultipleVariablesPostprocessor::initialize();
  _volume = 0;
}

void
ElementAverageMultipleValues::execute()
{
  ElementIntegralMultipleVariablesPostprocessor::execute();

  _volume += _current_elem_volume;
}

Real
ElementAverageMultipleValues::getValue()
{
  Real integral = ElementIntegralMultipleVariablesPostprocessor::getValue();

  gatherSum(_volume);

  return integral / _volume;
}

void
ElementAverageMultipleValues::threadJoin(const UserObject & y)
{
  ElementIntegralMultipleVariablesPostprocessor::threadJoin(y);
  const ElementAverageMultipleValues & pps = static_cast<const ElementAverageMultipleValues &>(y);
  _volume += pps._volume;
}
