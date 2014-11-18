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

#ifndef ELEMENTAVERAGEMULTIPLEVALUES_H
#define ELEMENTAVERAGEMULTIPLEVALUES_H

#include "ElementIntegralMultipleVariablesPostprocessor.h"
#include "EquationOfState.h"

//Forward Declarations
class ElementAverageMultipleValues;

template<>
InputParameters validParams<ElementAverageMultipleValues>();

/**
 * This postprocessor computes a volume integral of the specified variable.
 *
 * Note that specializations of this integral are possible by deriving from this
 * class and overriding computeQpIntegral().
 */
class ElementAverageMultipleValues : public ElementIntegralMultipleVariablesPostprocessor
{
public:
  ElementAverageMultipleValues(const std::string & name, InputParameters parameters);

  virtual void initialize();
  virtual void execute();
  virtual Real getValue();
  virtual void threadJoin(const UserObject & y);

protected:
    // Cell volume:
    Real _volume;
};

#endif
