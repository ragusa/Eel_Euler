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

#ifndef INTERNALENERGYAUX_H
#define INTERNALENERGYAUX_H

#include "AuxKernel.h"

class InternalEnergyAux;

template<>
InputParameters validParams<InternalEnergyAux>();

class InternalEnergyAux : public AuxKernel
{
public:

  InternalEnergyAux(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();

  bool _implicit;
  VariableValue & _rhoA;
  VariableValue & _rhouA_x;
  VariableValue & _rhouA_y;
  VariableValue & _rhouA_z;
  VariableValue & _rhoEA;
VariableValue & _area;
};

#endif //INTERNALENERGYAUX_H
