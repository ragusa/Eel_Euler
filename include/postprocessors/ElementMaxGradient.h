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

#ifndef ELEMENTMAXGRADIENT_H
#define ELEMENTMAXGRADIENT_H

#include "ElementPostprocessor.h"

//Forward Declarations
class ElementMaxGradient;

template<>
InputParameters validParams<ElementMaxGradient>();

/**
 * This postprocessor computes a volume integral of the specified variable.
 *
 * Note that specializations of this integral are possible by deriving from this
 * class and overriding computeQpIntegral().
 */
class ElementMaxGradient : public ElementPostprocessor
{
public:
  ElementMaxGradient(const std::string & name, InputParameters parameters);

  virtual void initialize();
  virtual void execute();
  virtual void finalize();
  virtual Real getValue();
  virtual void threadJoin(const UserObject & y);

protected:
    // Aux variable
    VariableGradient & _grad_press;
    // Variable:
    Real _value;
};

#endif
