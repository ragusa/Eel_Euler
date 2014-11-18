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

#ifndef NODALMINMULTIPLEVALUES_H
#define NODALMINMULTIPLEVALUES_H

#include "NodalVariablePostprocessor.h"
#include "EquationOfState.h"

class MooseVariable;

//Forward Declarations
class NodalMinMultipleValues;

template<>
InputParameters validParams<NodalMinMultipleValues>();

class NodalMinMultipleValues : public NodalVariablePostprocessor
{
public:
  NodalMinMultipleValues(const std::string & name, InputParameters parameters);

  virtual void initialize();
  virtual void execute();
  virtual Real getValue();
  virtual void threadJoin(const UserObject & y);

protected:
    Real _value;
    
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
