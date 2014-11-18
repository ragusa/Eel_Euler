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

#ifndef ELEMENTL1ERROR_H
#define ELEMENTL1ERROR_H

#include "ElementIntegralVariablePostprocessor.h"
// include "FunctionInterface.h"

class Function;

//Forward Declarations
class ElementL1Error;

template<>
InputParameters validParams<ElementL1Error>();

class ElementL1Error :
  public ElementIntegralVariablePostprocessor
  // public FunctionInterface
{
public:
  ElementL1Error(const std::string & name, InputParameters parameters);

  /**
   * Get the L2 Error.
   */
  virtual Real getValue();

protected:
  virtual Real computeQpIntegral();

  Function & _func;
};

#endif //ElementL1Error_H
