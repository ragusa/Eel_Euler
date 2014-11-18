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

#ifndef JUMPGRADIENTINTERFACE_H
#define JUMPGRADIENTINTERFACE_H

#include "InternalSideUserObject.h"

class JumpGradientInterface;

template<>
InputParameters validParams<JumpGradientInterface>();

/**
 *
 */
class JumpGradientInterface : public InternalSideUserObject
{
public:
  JumpGradientInterface(const std::string & name, InputParameters parameters);
  virtual ~JumpGradientInterface();

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
    VariableGradient & _grad_u;
    VariableGradient & _grad_u_neighbor;
    // Name of the variable storing the jump:
    std::string _jump_name;
    // Temporary variable:
    Real _value;
};

#endif /* JUMPGRADIENTINTERFACE_H */
