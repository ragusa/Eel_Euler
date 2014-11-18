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

#ifndef ELEMENTINTEGRALMULTIPLEVARIABLESPOSTPROCESSOR_H
#define ELEMENTINTEGRALMULTIPLEVARIABLESPOSTPROCESSOR_H

#include "ElementIntegralPostprocessor.h"
#include "MooseVariableInterface.h"
#include "EquationOfState.h"

//Forward Declarations
class ElementIntegralMultipleVariablesPostprocessor;

template<>
InputParameters validParams<ElementIntegralMultipleVariablesPostprocessor>();

/**
 * This postprocessor computes a volume integral of the specified variable.
 *
 * Note that specializations of this integral are possible by deriving from this
 * class and overriding computeQpIntegral().
 */
class ElementIntegralMultipleVariablesPostprocessor :
  public ElementIntegralPostprocessor,
  public MooseVariableInterface
{
public:
  ElementIntegralMultipleVariablesPostprocessor(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeQpIntegral();

    // Output type
    enum OutputType
    {
        RHOVEL2 = 0,
        RHOCVEL = 1,
        RHOC2 = 2
    };
    std::string _output_name;
    MooseEnum _output_type;
    // Variable
    MooseVariable & _var;
    // Conservative variables:
    VariableValue & _rhoA;
    VariableValue & _rhouA_x;
    VariableValue & _rhouA_y;
    VariableValue & _rhoEA;
    // Primitive variable:
    VariableValue & _area;
    // Equation of state:
    const EquationOfState & _eos;
};

#endif
