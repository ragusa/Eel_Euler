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

#ifndef PRESSUREAUX_H
#define PRESSUREAUX_H

#include "AuxKernel.h"
#include "EquationOfState.h"

//Forward Declarations
class PressureAux;

template<>
InputParameters validParams<PressureAux>();

class PressureAux : public AuxKernel
{
public:

  PressureAux(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();

  VariableValue & _rhoA;
  VariableValue & _rhouA_x;
  VariableValue & _rhouA_y;
  VariableValue & _rhouA_z;
  VariableValue & _rhoEA;
  VariableValue & _area;
  const EquationOfState & _eos;
};

#endif //PRESSUREAUX_H
