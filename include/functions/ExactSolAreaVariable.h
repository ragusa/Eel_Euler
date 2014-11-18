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

#ifndef EXACTSOLAREAVARIABLE_H
#define EXACTSOLAREAVARIABLE_H

#include "Function.h"
#include "FunctionInterface.h"
#include "EquationOfState.h"

class ExactSolAreaVariable;

template<>
InputParameters validParams<ExactSolAreaVariable>();

class ExactSolAreaVariable : public Function
{
public:
  ExactSolAreaVariable(const std::string & name, InputParameters parameters);

  virtual Real value(Real t, const Point & p);

  virtual RealVectorValue gradient(Real t, const Point & p);

protected:
    enum VariableType
    {
        DENSITY = 0,
        VELOCITY = 1,
        PRESSURE = 2
    };
    
    // Name of the variable the function will output.
    std::string _var_name;
    
    MooseEnum _var_type;
    
    // Stagnation variables:
    Real _Po;
    Real _To;
    
    // Back pressure
    Real _Pout;
    
    // Length of the computational domain
    Real _length;
    
    // Equation of state:
//    const EquationOfState & _eos;
};

#endif //EXACTSOLAREAVARIABLE_H
