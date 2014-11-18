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

#ifndef SMOOTHFUNCTION_H
#define SMOOTHFUNCTION_H

#include "InternalSideUserObject.h"

class SmoothFunction;

template<>
InputParameters validParams<SmoothFunction>();

/**
 *
 */
class SmoothFunction : public InternalSideUserObject
{
public:
  SmoothFunction(const std::string & name, InputParameters parameters);
  virtual ~SmoothFunction();

  virtual void initialize();
  virtual void execute();
  virtual void destroy();
  virtual void finalize();
  virtual void threadJoin(const UserObject & uo);

  Real getValue() const { return _value; }

protected:
    // Auxiliary system variable:
    AuxiliarySystem & _aux;
    // Gradient value:
    VariableValue & _u;
    VariableValue & _u_neighbor;
    // Name of the variable storing the jump:
    std::string _var_name;
    // Temporary variable:
    Real _value;
};

#endif /* SMOOTHFUNCTION_H */
